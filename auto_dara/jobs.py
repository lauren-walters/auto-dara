"""Runs searches, manual refinements and batch job files. Shared by the command line and the Python API."""

import csv
import dataclasses
import json
import re
import shutil
from pathlib import Path

from . import refine
from .cifs import build_pool
from .config import Settings

PATTERN_EXTS = {".xrdml", ".xy", ".raw", ".rasx"}
DICT_FIELDS = ("phase_params", "refinement_params", "search_phase_params", "search_refinement_params")
# where `-P PHASE:key=value` goes: manual refinement or the search
PER_PHASE_FIELDS = {"phase_params": "phase_overrides", "search_phase_params": "search_phase_overrides"}
_OPTIONAL_TYPES = {"rpb_threshold": float, "cpus": int, "custom_cifs": Path, "icsd_path": Path}


def for_search(overrides: dict) -> dict:
    """On a search, -P/-R (and those batch columns) set the search's parameters."""
    out = dict(overrides)
    for key in ("phase_params", "refinement_params"):
        if out.get(key):
            out["search_" + key] = out.pop(key)
    return out


def _value(text):
    if re.fullmatch(r"-?\d+", text):
        return int(text)
    if re.fullmatch(r"-?(\d+\.?\d*|\.\d+)([eE]-?\d+)?", text):
        return float(text)
    return text


def parse_params(items):
    """'k1=0_0^0.01' -> ({'k1': '0_0^0.01'}, {}); 'ZnO_186:k2=fixed' -> ({}, {'ZnO_186': {'k2': 'fixed'}})."""
    if isinstance(items, str):
        items = items.replace(";", " ").split()
    common, per_phase = {}, {}
    for item in items or []:
        key, sep, value = item.partition("=")
        if not sep:
            raise ValueError(f"parameter '{item}' should look like key=value")
        phase, _, key = key.strip().rpartition(":")
        (per_phase.setdefault(phase, {}) if phase else common)[key.strip()] = _value(value.strip())
    return common, per_phase


def settings_from(base: Settings, overrides: dict) -> Settings:
    s = dataclasses.replace(
        base,
        **{k: dict(getattr(base, k)) for k in DICT_FIELDS},
        phase_overrides={k: dict(v) for k, v in base.phase_overrides.items()},
        search_phase_overrides={k: dict(v) for k, v in base.search_phase_overrides.items()},
    )
    for f in dataclasses.fields(Settings):
        v = overrides.get(f.name)
        if v is None or v == "" or v == []:
            continue
        if f.name in PER_PHASE_FIELDS.values() and isinstance(v, dict):
            for phase, params in v.items():
                getattr(s, f.name).setdefault(phase, {}).update(params)
            continue
        if f.name in DICT_FIELDS:
            common, per_phase = (v, {}) if isinstance(v, dict) else parse_params(v)
            getattr(s, f.name).update(common)
            if per_phase and f.name not in PER_PHASE_FIELDS:
                raise ValueError(f"per-phase values ({', '.join(per_phase)}:...) only work for phase parameters (-P)")
            for phase, params in per_phase.items():
                getattr(s, PER_PHASE_FIELDS[f.name]).setdefault(phase, {}).update(params)
            continue
        default = getattr(Settings(), f.name)
        if isinstance(v, str):
            if f.name == "max_e_hull":
                v = None if v.lower() == "none" else float(v)
            elif isinstance(default, bool):
                v = v.strip().lower() in ("1", "true", "yes", "y")
            elif isinstance(default, tuple):
                v = tuple(v.replace(",", " ").split())
            elif isinstance(default, (int, float, Path)):
                v = type(default)(v)
            elif f.name in _OPTIONAL_TYPES:
                v = _OPTIONAL_TYPES[f.name](v)
        setattr(s, f.name, v)
    s.icsd_path = Path(s.icsd_path).expanduser() if s.icsd_path else None
    s.custom_cifs = Path(s.custom_cifs).expanduser() if s.custom_cifs else None
    s.databases = tuple(s.databases)
    return s


def load_settings(*folders) -> dict:
    for folder in folders:
        path = Path(folder) / "settings.json"
        if path.exists():
            return json.loads(path.read_text())
    return {}


def save_settings(s: Settings, folder: Path):
    folder.mkdir(parents=True, exist_ok=True)
    data = {k: (str(v) if isinstance(v, Path) else v) for k, v in dataclasses.asdict(s).items()}
    (folder / "settings.json").write_text(json.dumps(data, indent=2))


def find_patterns(items) -> list:
    out = []
    for item in map(Path, items):
        if item.is_dir():
            out += sorted(p for p in item.iterdir() if p.suffix.lower() in PATTERN_EXTS)
        elif item.exists():
            out.append(item)
        else:
            raise FileNotFoundError(f"pattern not found: {item}")
    return [p.resolve() for p in out]


def check_instrument(patterns, s):
    for p in patterns:
        if p.suffix.lower() == ".xrdml":
            head = p.read_text(errors="replace")[:20000].lower()
            if "aeris" in head and "aeris" not in s.instrument.lower():
                print(
                    f"WARNING: {p.name} was measured on an Aeris, but the instrument profile is {s.instrument}. "
                    "Use -i Aeris-fds-Pixcel1d-Medipix3 unless you know otherwise."
                )


def resolve_phase(name, *folders) -> Path:
    p = Path(name).expanduser()
    if p.suffix == ".cif" and p.exists():
        return p.resolve()
    stem = p.stem if p.suffix == ".cif" else name
    for folder in folders:
        lookup = refine.cif_lookup(Path(folder))
        for test in (lambda k: k == stem, lambda k: k.startswith(stem + "_"), lambda k: stem in k):
            hits = sorted(k for k in lookup if test(k))
            if len(hits) == 1:
                return lookup[hits[0]].resolve()
            if len(hits) > 1:
                raise ValueError(f"'{name}' matches several CIFs, be more specific:\n  " + "\n  ".join(hits))
    raise ValueError(f"no CIF matching '{name}'. Run a search first, or give a path to a .cif file.")


def _standard_name(cif: Path) -> Path:
    """Hand-added CIFs get the <formula>_<spacegroup>_custom_<name>.cif name the rest of auto-dara reads."""
    from pymatgen.core import Composition, Structure

    try:
        Composition(cif.stem.split("_")[0])
        return cif
    except Exception:
        pass
    st = Structure.from_file(str(cif))
    new = cif.with_name(f"{st.composition.reduced_formula}_{st.get_space_group_info()[1]}_custom_{cif.stem}.cif")
    cif.rename(new)
    print(f"  renamed {cif.name} -> {new.name}")
    return new


def run_search(patterns, s: Settings, out: Path) -> list:
    """Build the candidate pool once, then search each pattern. Each pattern gets its own folder in `out`."""
    check_instrument(patterns, s)
    first = out / refine.sample_name(patterns[0])
    if s.skip_filter:
        missing = [p.name for p in patterns if not any((out / refine.sample_name(p) / "cifs_considered").glob("*.cif"))]
        if missing:
            raise ValueError(
                f"--skip-filter needs an existing cifs_considered/ folder; none for {missing} in {out}. "
                "Run the search once without --skip-filter first."
            )
        pool = []
    else:
        pool = build_pool(s, first)
        if not pool:
            raise RuntimeError("no candidate structures found; widen the chemical system or add --cifs")
    rows = []
    for pattern in patterns:
        sample_dir = out / refine.sample_name(pattern)
        if s.skip_filter:
            sample_pool = [_standard_name(c) for c in sorted((sample_dir / "cifs_considered").glob("*.cif"))]
            print(f"\nReusing {len(sample_pool)} CIF(s) in {sample_dir / 'cifs_considered'} (filters skipped)")
        else:
            if sample_dir != first:
                shutil.rmtree(sample_dir / "cifs_considered", ignore_errors=True)
                shutil.copytree(first / "cifs_considered", sample_dir / "cifs_considered")
                shutil.copy2(first / "filter_report.md", sample_dir / "filter_report.md")
            sample_pool = [sample_dir / "cifs_considered" / p.name for p in pool]
        save_settings(s, sample_dir)
        try:
            rows += refine.search(pattern, sample_pool, s, out)
        except Exception as exc:
            print(f"  FAILED {pattern.name}: {type(exc).__name__}: {exc}")
        refine.write_summary(out)
        print(f"  saved in {sample_dir}")
    return rows


def run_manual(pattern: Path, phases, s: Settings, out: Path, add=()) -> dict:
    """Refine with the given phases; with none, start from the scan's cifs_best/. `add` appends phases."""
    check_instrument([pattern], s)
    sample_dir = out / refine.sample_name(pattern)
    if phases:
        cifs = [resolve_phase(n, sample_dir, out) for n in phases]
    else:
        cifs = sorted((sample_dir / "cifs_best").glob("*.cif"))
        if not cifs:
            raise ValueError(f"no phases given and no {sample_dir / 'cifs_best'} to start from; run a search first")
        print(f"Starting from cifs_best/: {', '.join(c.stem for c in cifs)}")
    cifs += [resolve_phase(n, sample_dir, out) for n in add]
    row = refine.manual(pattern, cifs, s, out)
    refine.write_summary(out)
    print(f"  saved in {sample_dir}")
    return row


def load_jobs(path: Path) -> list:
    path = Path(path)
    if path.suffix.lower() == ".json":
        jobs = json.loads(path.read_text())
    else:
        with path.open(newline="", encoding="utf-8-sig") as f:
            jobs = [{k.strip(): (v or "").strip() for k, v in row.items() if k} for row in csv.DictReader(f)]
    jobs = [j for j in jobs if str(j.get("pattern", "")).strip() and not str(j["pattern"]).startswith("#")]
    for j in jobs:
        j["pattern"] = str((path.parent / Path(j["pattern"]).expanduser()).resolve())
        if j.get("output"):
            j["output"] = str((path.parent / Path(j["output"]).expanduser()).resolve())
    return jobs


def _job_phases(job, column, base_dir):
    items = str(job.get(column) or "").replace(";", " ").split()
    return [str(base_dir / p) if p.endswith(".cif") and (base_dir / p).exists() else p for p in items]


def run_batch(path: Path, default_out: Path, base: Settings = None) -> list:
    """One row per pattern. Rows with `phases` (or `add`) are manual refinements, the rest are searches.
    `phases` = best uses the scan's cifs_best/.
    Search rows with identical settings and output folder share one candidate pool."""
    base = base or Settings()
    jobs = load_jobs(path)
    print(f"{len(jobs)} job(s) in {path}")
    status = []

    groups = {}
    for j in jobs:
        if j.get("phases") or j.get("add"):
            continue
        out = Path(j.get("output") or default_out)
        s = settings_from(base, for_search(j))
        key = (str(out), json.dumps(dataclasses.asdict(s), default=str, sort_keys=True))
        groups.setdefault(key, (out, s, []))[2].append(Path(j["pattern"]))

    for out, s, patterns in groups.values():
        print(f"\n##### search: {', '.join(p.name for p in patterns)}")
        try:
            rows = run_search(patterns, s, out)
            done = {r["sample"] for r in rows}
            status += [
                (
                    refine.sample_name(p),
                    "search",
                    "done" if refine.sample_name(p) in done else "no solution",
                    out / refine.sample_name(p),
                )
                for p in patterns
            ]
        except Exception as exc:
            print(f"  FAILED: {type(exc).__name__}: {exc}")
            status += [
                (refine.sample_name(p), "search", f"failed: {exc}", out / refine.sample_name(p)) for p in patterns
            ]

    for j in jobs:
        if not (j.get("phases") or j.get("add")):
            continue
        pattern = Path(j["pattern"])
        out = Path(j.get("output") or default_out)
        sample = refine.sample_name(pattern)
        print(f"\n##### manual: {pattern.name}")
        try:
            s = settings_from(settings_from(base, load_settings(out / sample)), j)
            phases = [p for p in _job_phases(j, "phases", path.parent) if p.lower() != "best"]
            run_manual(pattern, phases, s, out, add=_job_phases(j, "add", path.parent))
            status.append((sample, "manual", "done", out / sample))
        except Exception as exc:
            print(f"  FAILED: {type(exc).__name__}: {exc}")
            status.append((sample, "manual", f"failed: {exc}", out / sample))

    print("\n##### batch finished")
    width = max([len(r[0]) for r in status] + [6])
    for sample, kind, state, folder in status:
        print(f"  {sample:<{width}}  {kind:<12} {state:<12} {folder}")
    for out in sorted({str(r[3].parent) for r in status}):
        print(f"  summary: {Path(out) / 'summary.csv'}")
    return status
