"""Runs dara's phase search and manual refinements, ranks them by FOM, and saves the outputs."""

import csv
import hashlib
import json
import logging
import os
import re
import shutil
import sys
import tempfile
import time
import warnings
from contextlib import contextmanager
from datetime import datetime
from importlib import metadata
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
    "json",
]

# CIF names are {formula}_{space group number}_{icsd|cod|custom}_{id}, see cifs.py
_PHASE_NAME = re.compile(r"^(?P<formula>[^_]+)_(?P<sg>\d+)_(?P<db>icsd|cod|custom)_(?P<id>.+)$")
FOM_PARTS = {
    "gof": "gof_score",
    "missing_peaks": "missing_peaks_score",
    "extra_peaks": "extra_peaks_score",
    "intensity_mismatch": "intensity_mismatch_score",
    "element_presence": "element_presence_score",
    "lattice": "lattice_strain_score",
}


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


def _version(package):
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _number(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _value_and_uncertainty(x):
    """dara gives a refined number either as a plain value or as (value, uncertainty)."""
    if isinstance(x, (tuple, list)):
        return float(x[0]), (float(x[1]) if len(x) > 1 and x[1] is not None else None)
    return (None, None) if x is None else (float(x), None)


def _lattice(phase):
    """Lattice in angstrom and degrees. dara reports lengths in nm and leaves out values that the
    space group fixes (e.g. 90 degree angles, b = a for tetragonal); those are filled in here."""
    sg = phase.spacegroup_no or 0
    lat, unc = {}, {}
    for axis in ("a", "b", "c", "alpha", "beta", "gamma"):
        value, error = _value_and_uncertainty(getattr(phase, axis, None))
        scale = 10 if axis in "abc" else 1
        lat[axis] = None if value is None else round(value * scale, 5)
        unc[axis] = None if error is None else round(error * scale, 5)
    if lat["b"] is None and sg >= 75:
        lat["b"] = lat["a"]
    if lat["c"] is None and sg >= 195:
        lat["c"] = lat["a"]
    if 143 <= sg <= 194 and lat["alpha"] is not None:  # rhombohedral axes: all angles equal
        lat["beta"] = lat["beta"] if lat["beta"] is not None else lat["alpha"]
        lat["gamma"] = lat["gamma"] if lat["gamma"] is not None else lat["alpha"]
    for axis in ("alpha", "beta", "gamma"):
        if lat[axis] is None and sg >= 3:
            lat[axis] = 120.0 if (axis == "gamma" and 143 <= sg <= 194) else 90.0
    return lat, unc


def _phase_records(result):
    weights = result.get_phase_weights()
    phases = []
    for name, phase in result.lst_data.phases_results.items():
        parts = _PHASE_NAME.match(name)
        lattice, uncertainty = _lattice(phase)
        phases.append(
            {
                "name": name,
                "formula": parts["formula"] if parts else None,
                "spacegroup_number": phase.spacegroup_no,
                "spacegroup_symbol": phase.hermann_mauguin,
                "database": parts["db"] if parts else None,
                "database_id": parts["id"] if parts else None,
                "weight_fraction": round(float(weights[name]), 4),
                "lattice": lattice,
                "lattice_uncertainty": uncertainty,
            }
        )
    return sorted(phases, key=lambda p: -p["weight_fraction"])


def _fit_header(pattern, label, solution_found=True):
    return {
        "scan_file": Path(pattern).name,
        "scan_sha256": _sha256(pattern),
        "label": label,
        "kind": "manual" if label == "manual" else "automated",
        "solution_found": solution_found,
        "made_by": {"auto_dara": _version("auto-dara"), "dara": _version("dara-xrd")},
        "made_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    }


def _write_fit_json(prefix: Path, pattern, result, row: dict):
    """One small, readable JSON per fit with everything needed to report it without the .pkl."""
    lst = result.lst_data
    record = {
        **_fit_header(pattern, row["label"]),
        "rwp": round(lst.rwp, 3),
        "rexp": round(lst.rexp, 3),
        "goodness_of_fit": round(lst.rwp / lst.rexp, 3) if lst.rexp else None,
        "fom": _number(row.get("fom")),
        "fom_parts": {name: _number(row.get(col)) for col, name in FOM_PARTS.items()},
        "fom_note": "fom_parts are 0-1 scores (1 = best) that multiply to the FOM",
        "lattice_units": "angstrom and degrees",
        "phases": _phase_records(result),
        "files": {k: row.get(k) or None for k in ("png", "html", "pkl", "cifs")},
        "settings_file": "settings.json",
    }
    prefix.with_suffix(".json").write_text(json.dumps(record, indent=2))
    row["json"] = f"{prefix.name}.json"


def _save(result, fig, prefix: Path, cif_lookup: dict, row: dict, pattern, cif_dir_name=None):
    joblib.dump(result, prefix.with_suffix(".pkl"))
    fig.update_layout(title=f"{prefix.name}  Rwp {row['rwp']}  FOM {row['fom']}")
    # load plotly.js from the web instead of embedding it: ~1 MB per file instead of ~6 MB
    fig.write_html(prefix.with_suffix(".html"), include_plotlyjs="cdn")
    try:
        fig.write_image(prefix.with_suffix(".png"), width=1400, height=800)
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
    _write_fit_json(prefix, pattern, result, row)


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
        phases=pool,
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
        note = sample_dir / f"{sample}_no_solution.json"
        note.write_text(json.dumps(_fit_header(pattern, "rank0", solution_found=False), indent=2))
        _write_results(sample_dir, [{"sample": sample, "label": "rank0", "phases": "NO SOLUTION", "json": note.name}])
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
            pattern,
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
            if stem == prefix or stem.startswith(prefix + "_"):
                params[stem].update(extra)
    unused = [k for k in s.phase_overrides if not any(st == k or st.startswith(k + "_") for st in params)]
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
    _save(result, result.visualize(), sample_dir / f"{sample}_manual", lookup, row, pattern)
    _write_results(sample_dir, [row])
    print(f"  manual: FOM {row['fom'] or 'n/a'}  Rwp {row['rwp']}  {row['phases']}")
    return row
