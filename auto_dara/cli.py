"""auto-dara: ex situ XRD phase identification and refinement with dara.

auto-dara search scan.xrdml -p "La2O3, Fe2O3" -t LaFeO3 --temp "1100 C" -o results
auto-dara manual scan.xrdml LaFeO3_62 Fe2O3_167 -o results
auto-dara jobs.csv                      (many searches/refinements from one file)
auto-dara candidates -p "La2O3, Fe2O3" -o results
auto-dara pretty results/scan/scan_rank0.pkl         (publication-style PDF + PNG)
auto-dara instruments
"""

import argparse
import os
import sys
from pathlib import Path

from .config import Settings

COMMANDS = {"search", "manual", "batch", "candidates", "pretty", "instruments", "params", "-h", "--help"}


def _parsers():
    chem = argparse.ArgumentParser(add_help=False)
    g = chem.add_argument_group("sample")
    g.add_argument("-p", "--precursors", help='comma-separated formulas, e.g. "BaCO3, Fe2O3"')
    g.add_argument("-s", "--system", help="chemical system, e.g. Ba-Fe-O (added to the precursors' elements)")
    g.add_argument("-e", "--extra-elements", help='extra elements to search, e.g. "C, H" for air-exposed samples')
    g.add_argument("-t", "--target", help="expected product (context for the LLM)")
    g.add_argument(
        "-m",
        "--method",
        help='how it was made, e.g. "sol-gel"; "none" = unreacted precursors '
        "(default: solid-state synthesis in a box furnace)",
    )
    g.add_argument("--temp", dest="temperature", help='temperature and time, e.g. "1000 C for 12 h, furnace cooled"')
    g.add_argument("--atmosphere", help="synthesis atmosphere (default: air)")

    cand = argparse.ArgumentParser(add_help=False)
    g = cand.add_argument_group("candidate structures")
    g.add_argument("--db", dest="databases", nargs="+", choices=["icsd", "cod"], help="databases (default: icsd cod)")
    g.add_argument("--icsd", dest="icsd_path", type=Path, help="local ICSD folder of icsd_NNNNNN.cif files")
    g.add_argument("--cifs", dest="custom_cifs", type=Path, help="folder of your own CIFs to include")
    g.add_argument("--max-ehull", dest="max_e_hull", help="eV/atom; 'none' disables (default 0.1)")
    g.add_argument(
        "--no-llm",
        dest="use_llm",
        action="store_false",
        default=None,
        help="skip the LLM feasibility filter (use this if you have no Gemini API key)",
    )
    g.add_argument(
        "--skip-filter",
        dest="skip_filter",
        action="store_true",
        default=None,
        help="reuse each scan's existing cifs_considered/ folder (add your own CIFs there) and only rerun "
        "the phase search",
    )
    g.add_argument("--no-uniqueness", dest="use_uniqueness", action="store_false", default=None)

    ref = argparse.ArgumentParser(add_help=False)
    g = ref.add_argument_group("instrument")
    g.add_argument("-i", "--instrument", help="dara instrument profile (see `auto-dara instruments`)")
    g.add_argument("-w", "--wavelength", help="Cu, Co, Cr, Fe, Mo, or a number in angstrom (default Cu)")

    prm = argparse.ArgumentParser(add_help=False)
    g = prm.add_argument_group("refinement parameters (repeatable; see `auto-dara params`)")
    g.add_argument(
        "-P",
        "--param",
        dest="phase_params",
        action="append",
        metavar="KEY=VALUE",
        help="phase parameter, e.g. -P k1=0_0^0.01; prefix a phase to change only it: -P ZnO_186:k2=fixed",
    )
    g.add_argument(
        "-R",
        "--refine-param",
        dest="refinement_params",
        action="append",
        metavar="KEY=VALUE",
        help="whole-pattern parameter, e.g. -R eps2=0_-0.05^0.05",
    )

    srch = argparse.ArgumentParser(add_help=False)
    g = srch.add_argument_group("search")
    g.add_argument("--max-phases", type=int, help="most phases in one solution (default 5)")
    g.add_argument("--cases", dest="n_cases", type=int, help="solutions to save and rank (default 3)")
    g.add_argument("--rpb", dest="rpb_threshold", type=float, help="min Rpb gain to accept a phase (default: auto)")
    g.add_argument("--cpus", type=int, help="CPU cores for the search")

    fom = argparse.ArgumentParser(add_help=False)
    g = fom.add_argument_group("figure of merit")
    g.add_argument(
        "--no-element-check",
        dest="check_elements",
        action="store_false",
        default=None,
        help="do not require every precursor element to appear in a phase",
    )
    g.add_argument(
        "--optional-elements",
        help='elements that may appear in no phase, e.g. "K" for a K2CO3 flux or "Mg" for a dopant',
    )

    out = argparse.ArgumentParser(add_help=False)
    out.add_argument(
        "-o",
        "--out",
        type=Path,
        default=Path(os.environ.get("AUTO_DARA_RESULTS", "results")),
        help="output folder (default: $AUTO_DARA_RESULTS or ./results)",
    )

    root = argparse.ArgumentParser(
        prog="auto-dara", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = root.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser(
        "search", parents=[chem, cand, ref, prm, srch, fom, out], help="build candidates, search, rank by FOM"
    )
    p.add_argument("patterns", nargs="+", help="pattern files or folders (.xrdml .xy .raw .rasx)")
    p = sub.add_parser("candidates", parents=[chem, cand, out], help="only build and filter the candidate CIFs")
    p = sub.add_parser("manual", parents=[chem, ref, prm, fom, out], help="refine a pattern with phases you choose")
    p.add_argument("pattern", help="pattern file")
    p.add_argument(
        "phases",
        nargs="*",
        help="CIF paths, or names of CIFs from a previous search in --out. "
        "Leave out to use the scan's cifs_best/ (the best search solution)",
    )
    p.add_argument("--add", nargs="+", default=[], metavar="PHASE", help="phases to add, e.g. to test the target")
    p = sub.add_parser("batch", parents=[out], help="run every job in a .csv or .json file (see README)")
    p.add_argument("jobs", type=Path, help="job file; `auto-dara jobs.csv` is a shortcut for this")
    p = sub.add_parser("pretty", help="publication-style plot (PDF + high-res PNG) of saved refinements")
    p.add_argument("pkls", nargs="+", type=Path, help="result .pkl files, e.g. results/scan/scan_rank0.pkl")
    p.add_argument("--xlim", nargs=2, type=float, metavar=("MIN", "MAX"), help="2-theta range (default: all)")
    p.add_argument("-o", "--out", type=Path, default=Path("."), help="output folder (default: current folder)")
    sub.add_parser("instruments", help="list the instrument profiles you can pass to -i")
    sub.add_parser("params", help="show the default refinement parameters and what they mean")
    return root


def _instruments():
    from .instruments import list_profiles

    print("\n".join(list_profiles()))
    print("\nUse one with -i NAME. Default: Aeris-fds-Pixcel1d-Medipix3.")


def _pretty(args):
    from .pretty import save_pretty

    for pkl in args.pkls:
        try:
            pdf, png = save_pretty(pkl, args.out, args.xlim)
        except (FileNotFoundError, ValueError, ModuleNotFoundError) as exc:
            sys.exit(f"Error: {pkl}: {exc}")
        print(f"{pdf}\n{png}")


def _params():
    from . import config as c

    rows = [(k, c.DEFAULT_SEARCH_PHASE_PARAMS[k], v) for k, v in c.DEFAULT_PHASE_PARAMS.items()]
    rows += [(k, c.DEFAULT_SEARCH_REFINEMENT_PARAMS[k], v) for k, v in c.DEFAULT_REFINEMENT_PARAMS.items()]
    print(PARAMS_HELP.format(table="\n".join(f"  {k:<15}{str(a):<15}{b}" for k, a, b in rows)))


PARAMS_HELP = """Defaults. -P KEY=VALUE sets a phase parameter and -R KEY=VALUE sets eps1/eps2, for the command you run.
  parameter      search         manual
{table}

The search uses dara's own defaults, which its phase-acceptance test is tuned for. Manual refinement uses
slightly looser, tested values for a better final fit.

Values are BGMN refinement ranges, written start_min^max: "0_0^0.005" starts at 0 and may refine between
0 and 0.005. "fixed" holds a parameter at zero. lattice_range is the fraction the lattice may stretch (0.03 = 3%).

  k1, k2          strain-related peak broadening (see the BGMN manual)
  b1              crystallite-size broadening
  lattice_range   how far lattice parameters may move from the CIF
  gewicht         preferred orientation: SPHAR0 = none, SPHAR2/SPHAR4 = increasing texture freedom
  rp              BGMN peak-shape function
  eps1            zero-point shift of 2-theta
  eps2            sample displacement (height error)

Loosen with care: wide ranges let a wrong phase stretch and broaden until it fits anything.
"""


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] not in COMMANDS and Path(argv[0]).suffix.lower() in (".csv", ".json"):
        argv = ["batch"] + argv
    args = _parsers().parse_args(argv)
    if args.cmd == "instruments":
        return _instruments()
    if args.cmd == "params":
        return _params()
    if args.cmd == "pretty":
        return _pretty(args)
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass

    from . import jobs
    from .cifs import build_pool
    from .refine import sample_name

    out = args.out.resolve()
    try:
        if args.cmd in ("batch", "search", "manual"):
            from .rosetta import ensure_rosetta

            ensure_rosetta()
        if args.cmd == "batch":
            jobs.run_batch(args.jobs.resolve(), out)
        elif args.cmd == "search":
            patterns = jobs.find_patterns(args.patterns)
            base = Settings()
            if args.skip_filter:
                base = jobs.settings_from(base, jobs.load_settings(out / sample_name(patterns[0])))
            s = jobs.settings_from(base, jobs.for_search(vars(args)))
            jobs.run_search(patterns, s, out)
            print(f"\nDone. Each pattern has its own folder in {out}; summary.csv lists every solution.")
        elif args.cmd == "candidates":
            s = jobs.settings_from(Settings(), vars(args))
            jobs.save_settings(s, out / "candidates")
            build_pool(s, out / "candidates")
        elif args.cmd == "manual":
            pattern = jobs.find_patterns([args.pattern])[0]
            saved = jobs.load_settings(out / sample_name(pattern), out)
            s = jobs.settings_from(jobs.settings_from(Settings(), saved), vars(args))
            jobs.run_manual(pattern, args.phases, s, out, add=args.add)
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        sys.exit(f"Error: {exc}")


if __name__ == "__main__":
    main()
