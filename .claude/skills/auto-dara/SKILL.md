---
name: auto-dara
description: Run auto-dara to identify phases and Rietveld-refine ex situ powder XRD patterns (.xrdml, .xy, .raw, .rasx) with dara. Use when asked to analyse, refine or phase-identify XRD scans, run or check a phase search, do a manual refinement, test whether a phase is present, set up a batch job file, pick an instrument profile, or make a publication-quality refinement plot.
---

# auto-dara

Command-line tool in this repo. It finds candidate structures (local ICSD + COD), filters them (gases, E_hull,
structural duplicates, LLM chemical feasibility), runs dara's phase search, ranks solutions by figure of merit
(FOM), and saves everything. `README.md` is the user manual; `tutorial/README.md` is a worked example with known
answers.

## Setup

- Environment: conda env `auto-dara`, **Python 3.10** (`conda create -n auto-dara python=3.10 -y`, then
  `pip install -e .`). Don't use newer Python: package pins (`numpy<2`, `scipy<1.18`) matter.
- `.env` in the repo root (gitignored, never print or commit it): `GEMINI_API_KEY`, `ICSD_PATH`,
  `AUTO_DARA_RESULTS`. No `ICSD_PATH` means COD only. No key means use `--no-llm`.

## Commands

```bash
auto-dara search SCAN... -p "PRECURSORS" [-t TARGET] [-m METHOD] [--temp "T, time"] [--atmosphere air] -o OUT
auto-dara manual SCAN [PHASES...] [--add PHASE...] [-P KEY=VALUE] [-R KEY=VALUE] -o OUT
auto-dara JOBS.csv                     # batch: one row per scan; rows with `phases` are manual refinements
auto-dara candidates -p "..." -o OUT   # build and filter the pool only
auto-dara pretty PKL... [--xlim MIN MAX] [-o OUT]   # publication plot: <pkl>_pretty.pdf + 6000x3000 .png
auto-dara instruments                  # built-in instrument profiles for -i
auto-dara params                       # refinement parameter defaults and meanings
```

- `search --skip-filter` reruns only the phase search on the scan's existing `cifs_considered/` (use after
  adding or removing CIFs there). Settings are remembered from the first run.
- `manual` with no phases refines the scan's `cifs_best/`. `--add` puts extra phases on top. Phase names are
  CIF filename prefixes from the scan folder (e.g. `Sr2Nb2O7_36`); an ambiguous name lists the matches.
- `-P` changes phase parameters (`lattice_range`, `k1`, `k2`, `b1`, `gewicht`, `rp`), `-R` whole-pattern ones
  (`eps1`, `eps2`). A `PHASE:key=value` prefix works in manual refinement only.

## Output, per scan: `OUT/<scan>/`

`<scan>_rank0.*` (best by FOM; dara numbers from 0), `_rank1.*`, ... ; `<scan>_manual.*` (latest manual
refinement, replaced each time); `cifs_best/`; `cifs_considered/`; `filter_report.md`; `results.csv`;
`settings.json`. Plus `OUT/summary.csv` across scans. Each result is `.html` (interactive, loads plotly.js from the web), `.png`, `.json`
(phases, weight fractions, lattice in angstrom, Rwp, FOM) and `.pkl` (`joblib.load`). A search with no
solution writes `<scan>_no_solution.json`. Don't invent other output files or folders.

## Rules that matter

1. **Instrument profile must match the scan.** Read it from the file (`grep -o "Diffractometer system=[A-Z]*"`
   on `.xrdml`), never from another project's config. Default is `Aeris-fds-Pixcel1d-Medipix3`; README groups
   the built-in profiles by instrument. A wrong profile gives worse fits and much slower searches. Don't
   build custom profiles.
2. **Speed check (a diagnostic, not a limit).** auto-dara has no time limit; run time scales with the number of
   candidates. As a benchmark, ~30 candidates take 15 s to a few minutes. If a search that small runs past
   ~5 min, something is likely wrong (wrong instrument, CPU oversubscription): investigate rather than wait.
   Large pools (big systems, `--no-llm`) legitimately take much longer.
3. **Judge by FOM, not Rwp.** Rwp is background-dominated on noisy scans. FOM 0 means a required element is
   missing or the lattice drifted: open the plot.
4. **Search and manual use different defaults on purpose.** The search keeps dara's own parameters (looser ones
   let the main phase swallow minor phases); manual refinement uses slightly looser ones for the final fit.
   Don't unify them.
5. **Is phase X present?** `manual --add X`: real phases take weight and raise the FOM; absent ones refine to
   about 0 wt% with no FOM change.
6. **Ray:** never run `ray stop --force` (machine-wide). Stop a run with Ctrl-C or by killing its own process
   tree.
7. **Never modify dara or other shared checkouts.** Work around dara behaviour in this package.
8. Synthesis method and temperature drive the LLM's oxidation-state reasoning; pass them when known (`-m none`
   = unreacted precursors).

## Checking a run

Compare against the tutorial's known answers (run from `tutorial/`, outputs in `tutorial/results/`):
step 1 LaFeO3 FOM ~0.91; step 2 ZnAl2O4 + 8-10% ZnO; step 3 SrCO3 + Nb2O5 (target absent), manual FOM ~0.68.
