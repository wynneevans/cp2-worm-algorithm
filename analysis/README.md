# cp2-analysis

Analysis of the `M_<time>.dat` files written by the simulator. Python 3.11+, polars for data, seaborn for plots.

```zsh
# from the repo root (paths below are relative to it)
uv sync --project analysis
uv run --project analysis cp2-analysis <command> [FILE_OR_DIR ...] [options]
```
`FILE_OR_DIR` is an M file or a folder searched for finished runs (`M_<time>.dat` with its
`run_<time>.toml`; default `data`). Run parameters come from the run record; old flat `M_*.dat` files are
read through their long names when given as files. `--out` defaults to `<first input folder>/plots`; `--exact` to `tests/exact/exactCP2_results.csv`.
Plots and table rows are named by the run name and finish time (`raw_2x2_closed_20261028-114625.png`, ...).

## Commands
| Command | What it does | Writes |
| --- | --- | --- |
| `summary` | M and chi (charges 3 and 8) with jackknife errors; MC and exact drawn on the same axes | `m_summary.csv`, `susceptibility.png`, `magnetisation.png` |
| `raw` | every `--block` (100) measurements averaged into one point, per column | `raw_<name>_<time>.png` |
| `acf` | autocorrelation function and integrated autocorrelation time (Sokal window) | `acf_<name>_<time>.png`, `tau_int.csv` |
| `thermalisation` | chi3 and M3 against the number of rows skipped, and MSER-5 per column | `skipscan_<name>_<time>.png`, `mser.csv` |

Options: `--skip` (fraction of leading rows dropped, default 0.01; `summary`, `acf`), `--nblocks`
(jackknife blocks, default 50; `summary`, `thermalisation`), `--maxlag` (`acf`), `--block` (`raw`).

## Estimators
With `q3`, `q8` the temporal winding charges (doubled units: q3 = 2 Q3, q8 = (2/sqrt 3) Q8), `Lx` the
periodic length and `N = Lx*Ly`:

`M3 = <q3>/(2 Lx)`, `chi3 = beta <q3^2>/(4 N)`, `M8 = sqrt(3) <q8>/(2 Lx)`, `chi8 = 3 beta <q8^2>/(4 N)`.

The factors are constants in `observables.py`. They match the exact susceptibility in the three mu = 0 cases
in `tests/exact` (one beta, Lt = 4); check them again when other beta or Lt are analysed.

## Layout
```
src/cp2analysis/
  runs.py          RunInfo, file-name parsing, reading M_ files
  exact.py         exact results CSV
  stats.py         block_stats, acf, tau_int, mser (numpy only)
  observables.py   normalisations, summarise, skip_scan
  plotting.py      the only module touching matplotlib; theme, colours, save
  commands/        one module per subcommand
  cli.py           argparse entry point
tests/             pytest
```
Development: `uv run --project analysis pytest` and `... ruff check analysis`.
