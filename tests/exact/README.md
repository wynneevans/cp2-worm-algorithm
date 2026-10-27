# Exact results for Monte Carlo comparison

`exactCP2_results.csv` lists exact results for small lattices. The configs here run the
worm algorithm at the matching parameters so the two can be compared.

## The CSV
Columns: `lattice, boundary, beta_label, mu3, mu8, quantity, value, status, source, note`.
`quantity` is `susceptibility` or `magnetization`. Some rows appear twice (labelled file and
"plotted data"); they agree.

## Conventions used for the configs
- `beta_label` is the program's `beta`.
- `lattice` is `Lx x Ly`, with y the open direction for `y-open` (so `2x3` is `[2, 3, Lt]`).
- `Lt = 4`, so `dt = 4*beta/Lt = 0.05`. This is the only `dt` in `configs/probabilities.txt`,
  and the smallest `Lt` the program accepts. Smaller `dt` (for extrapolation towards the
  continuum) needs more rows in the probability table.

## mu3 = mu8 = 0
Three distinct cases, all at `beta = 0.05`:

| Config | Lattice | Boundary | Susceptibility | Magnetization |
| --- | --- | --- | --- | --- |
| `2x2_closed_mu0.toml` | 2x2 | closed | 0.0080440461 | 0 |
| `2x2_yopen_mu0.toml` | 2x2 | y-open | 0.0081171288 | 0 |
| `2x3_yopen_mu0.toml` | 2x3 | y-open | 0.008094515 | 0 |

Each runs 1,000,000 measurements (about 3 s) with a fixed seed, writing to `data/exact/`:
```zsh
build/release/src/wormCP2Mu2dAFMplaq tests/exact/2x2_closed_mu0.toml
build/release/src/wormCP2Mu2dAFMplaq tests/exact/2x2_closed_mu0.toml --seed 7   # independent repeat
```
The simulator does not thermalise, so discard the first measurements (10,000 is ample here)
when analysing. Each `M_*.dat` row is `Wx3, Wx8, Wy3, Wy8, Q3, Q8, l, lp, lm, tail`.

## Not done yet
There is no estimator that turns `M_*.dat` into a susceptibility or magnetization, and the
normalisation of the exact values (per site, per `beta`, ...) is not recorded here. As a
check, the mean of `Q3` and `Q8` should be consistent with 0 at `mu3 = mu8 = 0`, matching
the exact magnetization of 0. These runs are not part of `ctest`.
