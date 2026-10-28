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
