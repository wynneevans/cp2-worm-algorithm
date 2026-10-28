# cp2-analysis

Analysis of the `M_<time>.dat` files written by the simulator. Python 3.11+, polars for data, seaborn for plots.

```zsh
# from the repo root (paths below are relative to it)
uv sync --project analysis
uv run --project analysis cp2-analysis <command> [INPUT ...] [filters] [options]
```

## Choosing runs
An `INPUT` is the data folder, a parameter folder (`data/<bc>_<nmeas>_<beta>_<Lx>_<Ly>_<Lt>_mu3_<mu3>_mu8_<mu8>`)
or one `M_<time>.dat`; the default is `data`. Folders are searched for finished runs (`M_<time>.dat`
with its `run_<time>.toml`, where all parameters come from); unfinished `M_running_*` files are skipped.
Filters then narrow the selection by parameter. Each takes a value, a comma list or an inclusive range
`a:b` (`a:` and `:b` too): `--bc`, `--measurements`, `--beta`, `--mu3`, `--mu8`, `--lx`, `--ly`, `--lt`, `--lattice LxxLy[xLt]`.
A parameter without a filter can take any value.

```zsh
cp2-analysis summary data/y-open_100000_56.30_100_12_4504_mu3_0.010_mu8_0.000   # one parameter set
cp2-analysis summary --beta 56.3                       # the mu3 scan, plotted against mu3
cp2-analysis summary --beta 56.3 --mu3 0.01:0.03       # part of it
cp2-analysis summary --beta 56.3 --lx 100 --lt 4504    # every Ly and mu3: chi vs mu3, one series per Ly
cp2-analysis summary data/y-open_100000_56.30_*        # the same with a shell glob
cp2-analysis raw --beta 56.3 --mu3 0.02                # diagnostics of those runs
```
`--out` defaults to `<parameter folder>/plots` when one parameter set is selected, otherwise `plots/` in
the folder holding them (usually `data/plots`). `--exact` defaults to `tests/exact/exactCP2_results.csv`.

**Repeats.** A parameter folder can hold several runs. `summary` pools them by default: each run's first
`--skip` fraction is dropped, the rest of all runs joined and block-jackknifed (one result, `n_runs` in
the table). `--per-run` gives every run its own result instead, also skipping `--skip` of each.
`raw`, `acf` and `thermalisation` always look at single runs. Per-run files are named
`<run name>_<finish time>`, e.g. `raw_mu3_0.02_20261028-115713.png`.

**x axis.** `summary` plots against the first parameter that varies in the selection, in the order mu3,
mu8, beta, ly, lx, lt, measurements; each combination of the other varying parameters is its own series.
With nothing numeric varying it draws one tick per result. `--x {beta,mu3,mu8,lx,ly,lt,measurements,run}`
chooses. Exact values are drawn where the CSV has them.

## Commands
| Command | What it does | Writes |
| --- | --- | --- |
| `summary` | M and chi (charges 3 and 8) with jackknife errors, pooled per parameter set or `--per-run`; exact values on the same axes | `m_summary.csv`, `susceptibility.png`, `magnetisation.png` |
| `raw` | every `--block` (100) measurements averaged into one point, per column | `raw_<name>_<time>.png` |
| `acf` | autocorrelation function and integrated autocorrelation time (Sokal window) | `acf_<name>_<time>.png`, `tau_int.csv` |
| `thermalisation` | chi3 and M3 against the number of rows skipped, the same on `--chunks` (8) separate parts of the run, and MSER-5 per column | `skipscan_<name>_<time>.png`, `mser.csv`, `chunks.csv` |
| `check-exact` | pass/fail: M and chi (pooled as in `summary`) within `--max-pull` (4) errors of the exact values; `--typed` also uses the values typed into the notebook. Exit status 1 on a failure. Used by the `exact_*` ctest tests | none (prints a table) |

Options: `--skip` (fraction of each run's leading rows dropped, default 0.01; `summary`, `acf`), `--nblocks`
(jackknife blocks, default 50; `summary`, `thermalisation`), `--maxlag` (`acf`), `--block` (`raw`).

## Reading the thermalisation plots
- **Top row (skip scan)** is cumulative: each point uses every measurement after the skip, so neighbouring
  points share most of their data and a step between two points only reflects the slice removed between
  them. Compare steps with the error bars. A flat start means no thermalisation is needed.
- **Bottom row (chunks)** splits the run into separate parts with independent errors. If the run is
  stationary and the errors are right, the chunks agree: chi²/dof about their mean is near 1 (roughly
  0.3–2 for 8 chunks by chance). Well above that means a drift or, more often, an autocorrelation time
  longer than the run can measure, so the quoted errors are too small; repeats with other seeds tell
  which. `chunks.csv` holds the numbers and chi²/dof per run.

## Estimators
With `q3`, `q8` the temporal winding charges (doubled units: q3 = 2 Q3, q8 = (2/sqrt 3) Q8), `Lx` the
periodic length and `N = Lx*Ly`:

`M3 = <q3>/(2 Lx)`, `chi3 = beta <q3^2>/(4 N)`, `M8 = sqrt(3) <q8>/(2 Lx)`, `chi8 = 3 beta <q8^2>/(4 N)`.

The factors are constants in `observables.py`. They match the exact susceptibility in the three mu = 0 cases
in `tests/exact` (one beta, Lt = 4); check them again when other beta or Lt are analysed.

## Layout
```
src/cp2analysis/
  runs.py          RunInfo, finding finished runs, reading M_ files
  filters.py       --beta/--mu3/... selection
  exact.py         exact results CSV
  stats.py         block_stats, acf, tau_int, mser (numpy only)
  observables.py   normalisations, summarise, skip_scan
  plotting.py      the only module touching matplotlib; theme, colours, save
  commands/        one module per subcommand
  cli.py           argparse entry point
tests/             pytest (conftest.py builds a synthetic data folder)
```
Development: `uv run --project analysis pytest` and `... ruff check analysis`.
