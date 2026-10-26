"""Builds the filtered candidate CIF pool for a chemical system from ICSD, COD and your own CIFs."""

import itertools
import json
import os
import shutil
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .config import COD_CACHE, Settings

ELEMENTAL_GASES = {"H", "He", "N", "O", "F", "Ne", "Cl", "Ar", "Br", "Kr", "Xe", "Rn"}
COD_URL = "https://www.crystallography.net/cod/{}.cif"
COD_ATTEMPTS = 3
COD_RETRY_WAIT_S = 2  # times the attempt number
COD_WORKERS = 4  # COD drops connections under many parallel requests
SOURCE_PRIORITY = {"custom": 0, "icsd": 1, "cod": 2}


def source_of(name):
    return next((k for k in SOURCE_PRIORITY if f"_{k}_" in Path(name).name), "other")


def _common_gases():
    import dara

    return set(json.loads((Path(dara.__file__).parent / "data/common_gases.json").read_text()))


def _elements_of(formula):
    from pymatgen.core import Composition

    return {e.symbol for e in Composition(formula).elements}


def _norm(formula):
    from pymatgen.core import Composition

    try:
        return Composition(formula).reduced_formula.replace(" ", "")
    except Exception:
        return formula


def _db_entries(db, name, elements):
    info = db.preparsed_info
    out = []
    for r in range(1, len(elements) + 1):
        for sub in itertools.combinations(sorted(elements), r):
            for formula, code, sg, e_hull in info.get("-".join(sub), []):
                out.append({"source": name, "formula": _norm(formula), "code": str(code), "sg": sg, "e_hull": e_hull})
    return out


def _icsd_entries(s, elements):
    from dara.structure_db import ICSDDatabase

    from . import icsd

    e_hulls = {int(e["code"]): e["e_hull"] for e in _db_entries(ICSDDatabase(s.icsd_path), "icsd", elements)}
    out = []
    for e in icsd.entries_for_system(s.icsd_path, elements):
        code = e["icsd_id"]
        out.append(
            {
                "source": "icsd",
                "formula": icsd.format_reduced_formula(e["reduced"]),
                "code": code,
                "sg": e["sg"],
                "e_hull": e_hulls.get(int(code)) if code.isdigit() else None,
                "src": Path(s.icsd_path) / e["path"],
            }
        )
    return out


def _reject_reason(e, s, gases):
    els = _elements_of(e["formula"])
    if e["formula"] in gases or (len(els) == 1 and els <= ELEMENTAL_GASES):
        return "gas"
    if els <= {"C", "H", "O"}:
        return "C/H/O only"
    if "H" in els and "O" not in els:
        return "hydride"
    if s.max_e_hull is not None and e["e_hull"] is not None and e["e_hull"] > s.max_e_hull:
        return f"unstable (E_hull {e['e_hull'] * 1000:.0f} meV/atom)"
    return None


def _cod_reachable():
    import requests

    try:
        return requests.get(COD_URL.format(1000000), timeout=15).ok
    except Exception:
        return False


def _fetch_cod(code):
    """Returns (path, None) or (None, reason). Tries up to COD_ATTEMPTS times; a 404/410 is final."""
    import requests

    cached = COD_CACHE / f"{code}.cif"
    if cached.exists():
        return cached, None
    reason = None
    for attempt in range(1, COD_ATTEMPTS + 1):
        partial = cached.with_suffix(f".{os.getpid()}.part")  # other processes may share COD_CACHE
        try:
            r = requests.get(COD_URL.format(code), timeout=30)
            if r.status_code in (404, 410):
                return None, f"not on COD (HTTP {r.status_code})"
            r.raise_for_status()
            partial.write_bytes(r.content)
            partial.replace(cached)
            return cached, None
        except Exception as exc:
            partial.unlink(missing_ok=True)
            reason = f"download failed after {attempt} attempt(s): {type(exc).__name__}"
            if attempt < COD_ATTEMPTS:
                time.sleep(COD_RETRY_WAIT_S * attempt)
    return None, reason


def _download_cod(codes, problems):
    """Downloads `codes` into COD_CACHE; records why each one is missing in `problems`."""
    COD_CACHE.mkdir(parents=True, exist_ok=True)
    uncached = [c for c in codes if not (COD_CACHE / f"{c}.cif").exists()]
    if uncached and not _cod_reachable():
        problems.update({c: "COD unreachable" for c in uncached})
        codes = [c for c in codes if c not in uncached]
    with ThreadPoolExecutor(COD_WORKERS) as pool:
        fetched = dict(zip(codes, pool.map(_fetch_cod, codes), strict=True))
    problems.update({c: reason for c, (_, reason) in fetched.items() if reason})
    if problems:
        print(f"  warning: {len(problems)} COD structure(s) not available (listed in filter_report.md); "
              "a local COD copy avoids this, see `auto-dara mirror-cod`")
    return {c: path for c, (path, _) in fetched.items() if path}


def _stage_db(entries, s, stage_dir, log):
    from dara.structure_db import CODDatabase

    cod = CODDatabase()
    cod_codes = [e["code"] for e in entries if e["source"] == "cod"]
    cod_files, problems = {}, {}
    if cod_codes:
        if cod.local_copy_found:
            cod_files = {c: cod.get_file_path(c) for c in cod_codes if cod.get_file_path(c).exists()}
        # anything the local copy lacks (or all of it, without one) is downloaded
        missing = [c for c in cod_codes if c not in cod_files]
        if missing:
            cod_files.update(_download_cod(missing, problems))

    staged = []
    for e in entries:
        src = e.get("src") or cod_files.get(e["code"])
        if src is None or not Path(src).exists():
            reason = problems.get(e["code"], "file not available") if e["source"] == "cod" else "file not available"
            log.append(f"[missing] {e['formula']} {e['source']} {e['code']}: {reason}")
            continue
        dest = stage_dir / f"{e['formula']}_{e['sg']}_{e['source']}_{e['code']}.cif"
        shutil.copy2(src, dest)
        staged.append(dest)
    return staged


def _stage_custom(s, stage_dir, log):
    from pymatgen.core import Structure

    allowed = set(s.elements())
    staged = []
    for src in sorted(Path(s.custom_cifs).glob("*.cif")):
        try:
            st = Structure.from_file(str(src))
        except Exception as exc:
            log.append(f"[skip] custom {src.name}: unreadable ({exc})")
            continue
        if not {e.symbol for e in st.composition.elements} <= allowed:
            log.append(f"[skip] custom {src.name}: outside {s.chemsys()}")
            continue
        sg = st.get_space_group_info()[1]
        dest = stage_dir / f"{st.composition.reduced_formula}_{sg}_custom_{src.stem}.cif"
        shutil.copy2(src, dest)
        staged.append(dest)
    return staged


def _delete(paths):
    for p in paths:
        Path(p).unlink(missing_ok=True)


def build_pool(s: Settings, out_dir: Path) -> list:
    from dara.structure_db import CODDatabase

    elements = s.elements()
    if not elements:
        raise ValueError("give --precursors or --system so the chemical system is known")
    print(f"Chemical system: {'-'.join(elements)}")

    stage_dir = out_dir / "cifs_considered"
    shutil.rmtree(stage_dir, ignore_errors=True)
    stage_dir.mkdir(parents=True)
    log, report = [], [f"# Filter report: {'-'.join(elements)}", ""]

    entries = []
    if "icsd" in s.databases:
        if s.icsd_path and Path(s.icsd_path).exists():
            entries += _icsd_entries(s, elements)
        elif s.icsd_path:
            print(f"  warning: ICSD not found at {s.icsd_path}; using COD only")
        else:
            print("  no ICSD configured (set ICSD_PATH in .env or use --icsd); using COD only")
    if "cod" in s.databases:
        entries += _db_entries(CODDatabase(), "cod", elements)
    print(f"  {len(entries)} database entries")

    gases, best = _common_gases(), {}
    for e in entries:
        reason = _reject_reason(e, s, gases)
        if reason:
            log.append(f"[{reason}] {e['formula']} SG {e['sg']} {e['source']} {e['code']}")
            continue
        key = (e["formula"], e["sg"])
        rank = (SOURCE_PRIORITY[e["source"]], e["e_hull"] if e["e_hull"] is not None else 1.0)
        if key not in best or rank < best[key][0]:
            best[key] = (rank, e)
    kept = [v[1] for v in best.values()]
    print(f"  {len(kept)} after gas/hydride/E_hull filters and formula+space-group dedup")

    staged = _stage_db(kept, s, stage_dir, log)
    if s.custom_cifs:
        staged += _stage_custom(s, stage_dir, log)

    if s.use_uniqueness and staged:
        from .uniqueness import filter_unique

        res = filter_unique(staged, prefer=lambda p: SOURCE_PRIORITY.get(source_of(p), 3))
        _delete(res.discarded)
        report += ["## Structural duplicates", "", "```", res.report(), "```", ""]
        print(f"  uniqueness: {len(staged)} -> {len(res.kept)}")
        staged = sorted(res.kept)

    if s.use_llm and staged:
        try:
            from .feasibility import filter_feasible

            res = filter_feasible(staged, s)
            report += [res.report(s).replace("# Chemical feasibility", "## LLM chemical feasibility"), ""]
            if res.kept:
                _delete(res.discarded)
                print(f"  LLM feasibility: {len(staged)} -> {len(res.kept)}")
                staged = sorted(res.kept)
            else:
                print("  warning: LLM rejected everything; keeping full pool")
        except Exception as exc:
            print(f"  warning: LLM feasibility skipped ({exc})")

    report += ["## Dropped before staging (gases, hydrides, high E_hull, missing files)", "", "```", *log, "```"]
    (out_dir / "filter_report.md").write_text("\n".join(report))
    print(f"  {len(staged)} candidate CIF(s) in {stage_dir}")
    return staged
