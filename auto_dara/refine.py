"""Runs dara's phase search and manual refinements, ranks them by FOM, and saves the outputs."""

import csv
import logging
import os
import re
import shutil
import sys
import tempfile
import time
import warnings
from contextlib import contextmanager
from pathlib import Path

import joblib

from . import fom, peaks
from .config import Settings

_UUID = re.compile(r"_[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$")
FIELDS = [
    "sample",
    "label",
    "fom",
    "rwp",
    "phases",
    "gof",
    "missing_peaks",
    "extra_peaks",
    "intensity_mismatch",
    "element_presence",
    "lattice",
    "pkl",
    "html",
    "png",
    "cifs",
]


def sample_name(pattern):
    return _UUID.sub("", Path(pattern).stem)


def quiet():
    warnings.filterwarnings("ignore", category=UserWarning, module="dara")
    for name in ("dara", "dara.search.tree"):
        logging.getLogger(name).setLevel(logging.WARNING)
        logging.getLogger(name).propagate = False


@contextmanager
def _muted_stderr():
    sys.stderr.flush()
    saved = os.dup(2)
    with tempfile.TemporaryFile() as tmp:
        os.dup2(tmp.fileno(), 2)
        try:
            yield
        except BaseException:
            os.dup2(saved, 2)
            tmp.seek(0)
            sys.stderr.write(tmp.read().decode(errors="replace"))
            raise
        finally:
            os.dup2(saved, 2)
            os.close(saved)


def start_ray(cpus=None):
    import ray

    if ray.is_initialized():
        return
    with _muted_stderr():
        ray.init(
            num_cpus=cpus or int(os.environ.get("DARA_RAY_NUM_CPUS", 0)) or None,
            log_to_driver=False,
            include_dashboard=False,
            logging_level=logging.ERROR,
            runtime_env={"working_dir": None, "worker_process_setup_hook": peaks.install},
        )


def _save(result, fig, prefix: Path, cif_lookup: dict, row: dict, cif_dir_name=None):
    joblib.dump(result, prefix.with_name(f"{prefix.name}.pkl"))
    fig.update_layout(title=f"{prefix.name}  Rwp {row['rwp']}  FOM {row['fom']}")
    fig.write_html(prefix.with_name(f"{prefix.name}.html"))
    try:
        fig.write_image(prefix.with_name(f"{prefix.name}.png"), width=1400, height=800)
    except Exception as exc:
        print(f"  warning: PNG export failed ({exc})")
    row.update(pkl=f"{prefix.name}.pkl", html=f"{prefix.name}.html", png=f"{prefix.name}.png")
    if cif_dir_name:
        cif_dir = prefix.parent / cif_dir_name
        used = {n: Path(cif_lookup[n]).read_bytes() for n in result.get_phase_weights() if n in cif_lookup}
        shutil.rmtree(cif_dir, ignore_errors=True)
        cif_dir.mkdir()
        for name, data in used.items():
            (cif_dir / f"{name}.cif").write_bytes(data)
        row["cifs"] = cif_dir_name


def _row(sample, label, result, scores):
    weights = result.get_phase_weights()
    row = {
        "sample": sample,
        "label": label,
        "rwp": f"{result.lst_data.rwp:.3f}",
        "phases": ";".join(f"{n}:{w:.4f}" for n, w in sorted(weights.items(), key=lambda kv: -kv[1])),
    }
    row.update({k: "" if v is None else f"{v:.4f}" for k, v in scores.items()})
    return row


def _write_results(sample_dir: Path, rows):
    """Replace the rows being rewritten: all ranks after a search, the manual row after a manual refinement."""
    path = sample_dir / "results.csv"
    old = list(csv.DictReader(path.open())) if path.exists() else []
    new_ranks = any(r["label"].startswith("rank") for r in rows)
    labels = {r["label"] for r in rows}
    old = [r for r in old if r["label"] not in labels and not (new_ranks and r["label"].startswith("rank"))]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(old + rows, key=lambda r: (not r["label"].startswith("rank"), r["label"])))


def write_summary(out_dir: Path):
    rows = []
    for p in sorted(out_dir.glob("*/results.csv")):
        rows += list(csv.DictReader(p.open()))
    with (out_dir / "summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def cif_lookup(out_dir: Path) -> dict:
    """CIF name -> path under out_dir, preferring cifs_considered/ over the copies in cifs_best/."""
    files = sorted(out_dir.glob("**/*.cif"), key=lambda p: "cifs_considered" in p.parts)
    return {p.stem: p for p in files}


def matches_phase(stem: str, key: str) -> bool:
    """`key` names a phase by its CIF stem or a prefix of it (`Fe2O3_167`), or, for your own CIFs, by the file's
    original name (a custom CIF `my_phase.cif` is staged as `<formula>_<sg>_custom_my_phase.cif`)."""
    return stem == key or stem.startswith(key + "_") or stem.endswith("_custom_" + key)


def with_phase_overrides(pool: list, overrides: dict) -> list:
    """The pool for dara's search: a CIF that `overrides` names becomes a RefinementPhase carrying those
    parameters (on top of the search's phase parameters); every other CIF is passed as it is."""
    from dara.refine import RefinementPhase

    unused = [k for k in overrides if not any(matches_phase(Path(p).stem, k) for p in pool)]
    if unused:
        print(f"  warning: per-phase search parameters given for {unused}, which are not in the candidate pool")
    phases = []
    for p in pool:
        extra = {}
        for key, params in overrides.items():
            if matches_phase(Path(p).stem, key):
                extra.update(params)
        phases.append(RefinementPhase(path=Path(p), params=extra) if extra else p)
    return phases


def search(pattern: Path, pool: list, s: Settings, out_dir: Path) -> list:
    from dara.search.core import search_phases

    quiet()
    sample = sample_name(pattern)
    sample_dir = out_dir / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {sample}")
    print("  finding observed peaks (once per pattern, up to a few minutes on noisy scans)")
    n_peaks = fom.total_peaks(pattern, s)
    print(f"  searching {len(pool)} candidate(s); time grows with the number of candidates")
    start_ray(s.cpus)
    t0 = time.time()
    results = search_phases(
        pattern_path=pattern,
        phases=with_phase_overrides(pool, s.search_phase_overrides),
        wavelength=s.wavelength,
        instrument_profile=s.instrument,
        max_phases=s.max_phases,
        express_mode=True,
        enable_angular_cut=False,
        rpb_threshold=s.rpb_threshold,
        phase_params=dict(s.search_phase_params),
        refinement_params={**s.search_refinement_params, "n_threads": 1},
    )
    valid = [r for r in (results or []) if getattr(r, "refinement_result", None)]
    print(f"  search finished in {(time.time() - t0) / 60:.1f} min, {len(valid)} solution(s)")
    if not valid:
        print(f"  no solution found for {sample}")
        _write_results(sample_dir, [{"sample": sample, "label": "rank0", "phases": "NO SOLUTION"}])
        return []

    for old in [*sample_dir.glob(f"{sample}_rank*"), sample_dir / "cifs_best"]:
        shutil.rmtree(old) if old.is_dir() else old.unlink(missing_ok=True)
    lookup = cif_lookup(sample_dir)
    scored = []
    for case, r in enumerate(valid[: s.n_cases]):
        sc = fom.score(r.refinement_result, lookup, s.required_elements(), n_peaks, s)
        scored.append((case, r, sc))
    scored.sort(key=lambda t: (-(t[2].get("fom") if t[2].get("fom") is not None else -1), t[0]))

    rows = []
    for rank, (case, r, sc) in enumerate(scored):
        row = _row(sample, f"rank{rank}", r.refinement_result, sc)
        _save(
            r.refinement_result,
            r.visualize(),
            sample_dir / f"{sample}_rank{rank}",
            lookup,
            row,
            "cifs_best" if rank == 0 else None,
        )
        rows.append(row)
        print(f"  rank{rank} (dara case {case}): FOM {row['fom'] or 'n/a'}  Rwp {row['rwp']}  {row['phases']}")
    _write_results(sample_dir, rows)
    return rows


def manual(pattern: Path, cifs: list, s: Settings, out_dir: Path) -> dict:
    from dara.refine import RefinementPhase, do_refinement

    quiet()
    peaks.install()
    sample = sample_name(pattern)
    sample_dir = out_dir / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {sample}: manual refinement with {', '.join(Path(c).stem for c in cifs)}")

    params = {}
    for c in cifs:
        stem = Path(c).stem
        params[stem] = dict(s.phase_params)
        for prefix, extra in s.phase_overrides.items():
            if matches_phase(stem, prefix):
                params[stem].update(extra)
    unused = [k for k in s.phase_overrides if not any(matches_phase(st, k) for st in params)]
    if unused:
        raise ValueError(f"per-phase parameters given for {unused}, which are not in this refinement")

    work = sample_dir / "manual_work"
    shutil.rmtree(work, ignore_errors=True)
    result = do_refinement(
        pattern_path=pattern,
        phases=[RefinementPhase(path=c, params=params[Path(c).stem]) for c in cifs],
        wavelength=s.wavelength,
        instrument_profile=s.instrument,
        working_dir=work,
        refinement_params=dict(s.refinement_params),
    )
    shutil.rmtree(work, ignore_errors=True)

    lookup = {Path(c).stem: Path(c) for c in cifs}
    sc = fom.score(result, lookup, s.required_elements(), fom.total_peaks(pattern, s), s)
    row = _row(sample, "manual", result, sc)
    _save(result, result.visualize(), sample_dir / f"{sample}_manual", lookup, row)
    _write_results(sample_dir, [row])
    print(f"  manual: FOM {row['fom'] or 'n/a'}  Rwp {row['rwp']}  {row['phases']}")
    return row
