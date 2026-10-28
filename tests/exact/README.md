# Exact results for Monte Carlo comparison

`exactCP2_results.csv` lists exact results for small lattices. The configs here run the
worm algorithm at the matching parameters so the two can be compared.

## The CSV
Columns: `lattice, boundary, beta_label, mu3, mu8, quantity, value, status, source, note`.
`quantity` is `susceptibility` or `magnetization`. Some rows appear twice (labelled file and
"plotted data"); they agree.

## Conventions used for the runs
- `beta_label` is the program's `beta`.
- `lattice` is `Lx x Ly`, with y the open direction for `y-open` (so `2x3` is `[2, 3, Lt]`).
- `Lt = 4`, so `dt = 4*beta/Lt = 0.05`. This row is in `configs/probabilities.txt`, and `Lt = 4` is
  the smallest `Lt` the program accepts. Smaller `dt` (for extrapolation towards the continuum)
  needs more rows in the master probability table.

## mu3 = mu8 = 0
The three runs are in `configs/exact.toml`, all at `beta = 0.05`:

| Run name | Lattice | Boundary | Susceptibility | Magnetization |
| --- | --- | --- | --- | --- |
| `2x2_closed` | 2x2 | closed | 0.0080440461 | 0 |
| `2x2_yopen` | 2x2 | y-open | 0.0081171288 | 0 |
| `2x3_yopen` | 2x3 | y-open | 0.008094515 | 0 |

Each runs 1,000,000 measurements (about 3 s) with a fixed seed, writing to its parameter folder, e.g.
`data/closed_1000000_0.05_2_2_4_mu3_0.000_mu8_0.000/`:
```zsh
build/release/src/wormCP2Mu2dAFMplaq configs/exact.toml                       # all three
build/release/src/wormCP2Mu2dAFMplaq configs/exact.toml --only 2x2_closed --seed 7   # adds a repeat to the same folder
```
The simulator does not thermalise, so discard the first measurements (10,000 is ample here) when
analysing. Each `M_<time>.dat` row is `Wx3, Wx8, Wy3, Wy8, Q3, Q8, l, lp, lm, tail`.

Compare with the exact values using `cp2-analysis summary` (see the main README).

## Automatic checks in ctest
`exact_runs.toml` holds six fixed-seed runs (10^6 measurements each) that `ctest` compares with this CSV:
the three cases above, 2x3 y-open at mu3 = 0.03, 0.1 and 0.2 (exact values typed into the notebook) and
2x2 closed at mu8 = 0.1 and 0.2 (exact values from the plotted data).
`run_exact.cmake` runs one with the simulator, then `cp2-analysis check-exact <folder> --typed`, which
computes M and chi exactly as `summary` does and fails if |MC - exact| is more than 4 errors, or if the
CSV has no exact value. At mu3 = mu8 = 0 both charges are checked; otherwise the charge with the chemical
potential. The tests need uv and a synced analysis environment (`uv sync --project analysis`); if CMake
finds no `uv` they show as disabled. They take about 5 s each and are labelled `exact`: `ctest -L exact`
runs only them, `ctest -LE exact` skips them, `ctest -j4` runs them in parallel.

**What they can catch.** The errors are about 0.25% on chi and 0.0012 on M, so they catch broken updates,
wrong probability rows or normalisation slips, but not effects smaller than that: the exact chi barely
changes with mu3 here, and M at mu3 = 0.2 (0.0049) is only about 4 errors from zero.

**What is missing.** The master table only has dt = 0.05 rows with mu3 <= 0.4 (mu8 = 0) and mu8 <= 0.3
(mu3 = 0), so the exact results at mu3 or mu8 = 0.5 to 2 and at beta != 0.05 cannot be simulated. Adding those rows to
`configs/probabilities.txt`, a `[[run]]` to `exact_runs.toml` and its name to the list in
`tests/CMakeLists.txt` would add them.
