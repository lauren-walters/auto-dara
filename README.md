# auto-dara

Phase identification and Rietveld refinement of ex situ powder XRD patterns with
[dara](https://github.com/CederGroupHub/dara). Give it a pattern and the sample's chemistry; it finds candidate
structures, filters them, runs dara's phase search and ranks the solutions by figure of merit (FOM).

New to this? Work through `tutorial/README.md` after setup.

## Install

Use Python 3.10; newer versions pull in packages dara doesn't work with yet.

```bash
git clone https://github.com/lauren-walters/auto-dara.git
cd auto-dara
conda create -n auto-dara python=3.10 -y
conda activate auto-dara
pip install -e .
auto-dara --help
```

On an Apple Silicon Mac (M1 or later), dara's refinement engine (BGMN) is Intel-only and runs through Rosetta 2.
auto-dara checks for Rosetta when you run a search and offers to install it if it's missing. To install it yourself:
`softwareupdate --install-rosetta --agree-to-license`.

## Set up `.env`

From the `auto-dara` folder, run this with your own values filled in:

```bash
cat > .env << 'EOF'
GEMINI_API_KEY=paste-your-key-here
ICSD_PATH=/path/to/ICSD_2024/ICSD_2024_experimental_inorganic
AUTO_DARA_RESULTS=/path/to/results
EOF
```

Check that auto-dara reads it:

```bash
python -c "import auto_dara, os; [print(k, '=', os.environ.get(k) or 'not set') for k in ('ICSD_PATH', 'AUTO_DARA_RESULTS')]; print('GEMINI_API_KEY', 'set' if os.environ.get('GEMINI_API_KEY') else 'not set')"
```

Every setting is optional. To skip one, leave it empty (`ICSD_PATH=`).

| setting | what it is | if empty |
| --- | --- | --- |
| `GEMINI_API_KEY` | key for the LLM feasibility filter, from [Google AI Studio](https://aistudio.google.com/apikey) | add `--no-llm` to searches |
| `ICSD_PATH` | folder of ICSD `.cif` files (experimental only) | COD only |
| `AUTO_DARA_RESULTS` | default output folder | `./results` |

`.env` is read from the `auto-dara` folder or the folder you run from. `--icsd PATH` and `-o FOLDER` override it
for one run. The first search against a new ICSD folder takes a few minutes to index it.

## Search

```bash
auto-dara search scan.xrdml -p "La2O3, Fe2O3" -t LaFeO3 --temp "1100 C for 12 h" -o results
```

| option | meaning |
| --- | --- |
| `-p` / `-s` | precursors (`"La2O3, Fe2O3"`) or chemical system (`La-Fe-O`); at least one is required |
| `-e` | extra elements, e.g. `"C, H"` for carbonates or hydroxides |
| `-t`, `-m`, `--temp`, `--atmosphere` | target, synthesis method, temperature/time, atmosphere (default air); used by the LLM filter |
| `-m none` | no synthesis: you're measuring the precursors themselves |
| `-i` | instrument profile; must match your diffractometer (default `Aeris-fds-Pixcel1d-Medipix3`) |
| `--no-llm` | skip the LLM filter |
| `--skip-filter` | rerun only the phase search on the scan's existing `cifs_considered/` (add or delete CIFs there first) |
| `--optional-elements` | elements the FOM does not require in any phase, e.g. `"K"` for a K2CO3 flux or `"Mg"` for a dopant |
| `--no-element-check` | turn off the FOM's element check: no element is required in any phase |

Several patterns, or a folder, can be given at once if they share a chemical system. `auto-dara candidates ...`
builds the candidate pool without searching.

## Manual refinement

```bash
auto-dara manual scan.xrdml -o results                         # phases of the best search solution
auto-dara manual scan.xrdml --add Sr2Nb2O7_36 -o results       # best solution plus another phase
auto-dara manual scan.xrdml LaFeO3_62 Fe2O3_167 -o results     # exactly these phases
```

Phases are `.cif` paths or (partial) names of CIFs from an earlier search in the same `-o` folder. The result
is saved as `<scan>_manual.*`, replacing the previous one.

## Refinement parameters

`auto-dara params` lists the defaults and what they mean. The search uses dara's defaults; manual refinement uses
slightly looser ones. Change them with `-P` (phase parameters) and `-R` (`eps1`, `eps2`):

```bash
auto-dara manual scan.xrdml LaFeO3_62 Fe2O3_167 -o results -P k1=0_0^0.01 -P Fe2O3_167:lattice_range=0.05 -R eps2=0_-0.05^0.05
```

## Instrument profiles

`auto-dara instruments` lists them. Pick the one for your instrument model and slit setup (`fds` fixed
divergence slit, `ads` automatic).

| instrument | profiles |
| --- | --- |
| PANalytical Aeris | `Aeris-fds-Pixcel1d-Medipix3` |
| PANalytical X'Pert / Empyrean | `xpert-pixcel-0500`, `xpert-xcel-ads-10mm`, `xpert-xcel-ads-10mm-Ge`, `xpert-xcel-fds-0125`, `xpert-xcel-fds-0250`, `xpert-xcel-fds-1000`, `xpert-xcel-htk-ads-10`, `xpert-xcel-htk-fds-0125`, `PW3040-FDS-ADS-Xcelerator` |
| PANalytical Cubix | `cubix-ads-10mm`, `cubix-ads-15mm` |
| Philips PW1800 | `pw1800-ads-10mm`, `pw1800-fds` |
| Bruker D8 | `d8-fds-02-LynxEyeXE`, `d8-lynxeye-ads-1mm`, `d8-lynxeye-fds-02`, `d8-lynxeye-fds-05mm`, `d8-lynxeye-fds-06mm`, `d8-solxe-fds-0600`, `d8-solxe-vds-12mm`, `D8_6Div_4SS`, `LBL-d8-LynxEyeXE`, `RMS-D8-ADS-15-Glass-LynxEyeXE`, `RMS-D8-ADS-15-LynxEyeXE`, `RMS-D8-Capillary-500um-LynxEyeXE`, `RMS-D8-FDS-03-LynxEyeXE` |
| Bruker D2 | `d2-ssd160-fds-1` |
| Bruker X8 Apex2 | `x8-apex2-fds-10` |
| Siemens D5000 | `siemens-d5000-fds1mm`, `siemens-d5000-fds1mm-2`, `siemens-d5000-fds2mm` |
| Rigaku MiniFlex | `Rigaku-Miniflex`, `Rigaku-Miniflex-600-DTEXultra2-fds`, `Rigaku-Miniflex-gen5-Dtex-var`, `Rigaku-Miniflex-gen5-Dtex-varfix` |
| Rigaku SmartLab / Ultima | `Rigaku-SmartLab-CBO-BB-FDS05deg`, `rigaku-ultima` |
| Anton Paar XRDynamic 500 | `ucb-xrdynamic500-anton-paar` |
| Synchrotron | `synchrotron` |

## Publication plots

```bash
auto-dara pretty results/scan/scan_rank0.pkl                   # writes scan_rank0_pretty.pdf / .png in the current folder
auto-dara pretty results/*/*_rank0.pkl --xlim 10 60 -o figures
```

| option | meaning |
| --- | --- |
| `--xlim MIN MAX` | 2θ range to plot (default: the whole scan) |
| `-o` | output folder (default: the current folder) |

## Batch runs

```bash
auto-dara jobs.csv
```

One row per pattern. `pattern` is required (relative to the job file); other columns match the command-line
options (`precursors`, `system`, `target`, `temperature`, `instrument`, `output`, ...). Fill in `phases` (CIF
names, or `best`) to make a row a manual refinement. Empty cells use the defaults. See `tutorial/tutorial4.csv`.

## Results

```
results/
  summary.csv                   every solution for every scan
  <scan>/
    <scan>_rank0.html / .png    best solution (highest FOM); rank1, rank2 are next best
    <scan>_rank0.pkl            full dara result (joblib.load)
    <scan>_manual.*             latest manual refinement
    cifs_best/                  structures in the best solution
    cifs_considered/            every structure the search tried
    filter_report.md            what each filter removed and why
    results.csv                 FOM, Rwp, phases, weight fractions
    settings.json               exactly what was run
```

## Citing

auto-dara runs on [dara](https://github.com/CederGroupHub/dara), which drives the
[BGMN](http://www.bgmn.de) Rietveld program. If you use auto-dara in published work, please cite both:

```bibtex
@article{Fei2026,
  title     = {Dara: Automated Multiple-Hypothesis Phase Identification and Refinement from Powder X-ray Diffraction},
  author    = {Fei, Yuxing and McDermott, Matthew J. and Rom, Christopher L. and Wang, Shilong and Ceder, Gerbrand},
  journal   = {Chemistry of Materials},
  volume    = {38},
  number    = {3},
  pages     = {1364--1376},
  year      = {2026},
  month     = jan,
  publisher = {American Chemical Society (ACS)},
  issn      = {1520-5002},
  doi       = {10.1021/acs.chemmater.5c02820},
  url       = {http://dx.doi.org/10.1021/acs.chemmater.5c02820}
}

@article{BGMN,
  title     = {Handling of Unusual Instrumental Profiles by the BGMN Rietveld Program},
  author    = {Bergmann, J. and Friedel, P. and Kleeberg, R.},
  journal   = {Mater. Sci. Forum},
  volume    = {321--324},
  pages     = {192--197},
  year      = {2000},
  publisher = {Trans Tech Publications, Ltd.},
  issn      = {1662-9752},
  doi       = {10.4028/www.scientific.net/msf.321-324.192}
}
```

## License

MIT, see `LICENSE`.
