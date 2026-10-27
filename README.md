# CP(2) Worm Algorithm in C++

I wrote this code as a postdoc to perform the research carried out in this paper: https://arxiv.org/abs/1803.04767.

## Overview of the Worm Algorithm
We simulate a spin system which has two spatial dimensions and a time dimension. 

Hence a 3-dimensional lattice is used. On each lattice site there sits a CP(2) degree of freedom. The details of CP(2) - a mathematical symmetry group - are not important for those only interested in the software engineering aspects of this repository. It suffices to imagine a lattice of sites, with a (spin) value at each site that can take one of three values, let's say 0, 1 or 2.

We would like to update the current lattice configuration, a configuration with of 0s, 1s and 2s at each lattice site, so that we can generate thousands of different configurations. Once we have an ensemble of configurations we can then perform measurements one each configuration and by averaging the results we can learn about the physics of the CP(2) spin system.

To update the current configuration we use the worm algorithm:

i) First we randomly pick an initial time direction and starting site. At this point imagine the head and tail of the worm are both at the starting point.
ii) Then according to rules determined by the physics of the system, the value at the lattice site is updated and the head of the worm moves to another lattice site, this move can be:
a) spatially - if so the time direction of the worm also reverses
b) temporally
c) diagonally - spatially and temporally simultaneously.
iii) Step ii) is repeated until the head of the worm finds its tail, at which point the update is complete, and a new valid configuration has been produced.

## Repository layout
```
src/worm/          the worm algorithm simulator (wormCP2Mu2dAFMplaq)
src/tools/         analysis programs: analyseWorm, getWormXi, getXIutvNFiles
include/cp2worm/   shared headers (jackknife.h)
third_party/       ranlxd (Luescher's RANLUX generator, GPL) and toml++ (TOML parser, MIT)
tests/             CTest regression tests and golden output
tests/exact/       exact results and configs for Monte Carlo comparison
configs/           run configs (production.toml) and the probability table
scripts/           convert_probs.py, converts old probability files
analysis_with_python/  Python (uv, polars, seaborn) analysis of the output files
reference/         archive: original /dev/urandom-seeded simulator (not built) and the old analysis/ notebook
data/              simulation output (gitignored)
```

## Quick start
```zsh
# 1. Build (set Boost_ROOT to your Boost 1.74+ headers; see "Building" below)
cmake -S . -B build/release -DBoost_ROOT=/path/to/boost
cmake --build build/release

# 2. Check the build: runs the regression tests (a few seconds)
(cd build/release && ctest --output-on-failure)

# 3. Preview the production batch without running anything
build/release/src/wormCP2Mu2dAFMplaq configs/production.toml --dry-run

# 4. Run one small case, writing to data/test so nothing else is overwritten
build/release/src/wormCP2Mu2dAFMplaq configs/production.toml --only closed_b2.5_4x4x200 --output-dir data/test
```
Run it with no arguments, or with `--help`, to print the usage text.

## Building
Requires CMake 3.19+, a C++17 compiler and Boost 1.74+ (headers only).

```zsh
cmake -S . -B build/release -DBoost_ROOT=/path/to/boost
cmake --build build/release
(cd build/release && ctest --output-on-failure)
```

Machine-specific settings (Boost location, compiler flags) belong in a
`CMakeUserPresets.json`, which is gitignored. Copy the `local-base` pattern from
`CMakePresets.json`'s `release` preset and then run `cmake --preset=release-local`.

## Running a simulation
The simulator reads a TOML run file; `configs/production.toml` is the batch that used to be
`run_worm.zsh`:

```zsh
build/release/src/wormCP2Mu2dAFMplaq configs/production.toml            # all enabled runs
build/release/src/wormCP2Mu2dAFMplaq configs/production.toml --only closed_b5_8x8x400
build/release/src/wormCP2Mu2dAFMplaq configs/production.toml --dry-run   # show what would run
```

```toml
[defaults]
bc            = "closed"            # or "y-open"
measurements  = 100000
probabilities = "probabilities.txt" # relative to this file
output        = "../data"
mu3 = 0.0
mu8 = 0.0
sim = 1                             # tag in the output file names

[[run]]
name    = "b2.5_4x4x200"
beta    = 2.50
lattice = [4, 4, 200]               # Lx, Ly, Lt (Lt must be a multiple of 4)
# seed = 42                         # optional
# enabled = false                   # skip this run unless chosen with --only
```
### Command-line options
```
wormCP2Mu2dAFMplaq <run.toml> [options]
```
| Option | Effect |
| --- | --- |
| `--only NAME` | Run only the `[[run]]` called `NAME`, even if it has `enabled = false`. |
| `--output-dir DIR` | Write output to `DIR` instead of the file's `output`. |
| `--seed N` | Use seed `N` (1 to 2147483647) for every run. |
| `--dry-run` | Print the resolved runs and the probability row each will use, then exit. |
| `-h`, `--help` | Show the usage text. |

The program exits with status 1 and an `error:` message on bad input (config, probability
table or options), before any simulation starts.

Any `[defaults]` key can be overridden inside a `[[run]]`. Unknown keys are errors, so a
typo cannot silently fall back to a default. Everything is checked before the first run
starts.

**Seeds.** With no `seed` (and no `--seed N`) each run draws a random one, prints it, and
writes `run_<tag>.toml` next to its `.dat` files. Running that file reproduces the output
byte for byte: `wormCP2Mu2dAFMplaq <output>/run_<tag>.toml --output-dir elsewhere`.

**Probability table.** `configs/probabilities.txt` has a header line and one row per
`(dt, mu3, mu8)`, where `dt = 4*beta/Lt`:
```
# comments start with '#'
dt mu3 mu8 p01 p02 ... p48
0.05 0 0 1. 0 0.0253044 ...
```
Each run picks the row matching its own `4*beta/Lt`, `mu3` and `mu8`, and stops with the list
of available rows if there is none. To simulate a new `dt`, add a row for it. Older
headerless files (51 numbers per line) convert with
`scripts/convert_probs.py old.txt new.txt`.

### Output
Each run writes, into its output directory:

| File | Contents |
| --- | --- |
| `C_wormCP2AFMmu2dplaq_<tag>.dat` | Correlator measurements, comma-separated, read by `analysis/` and the tools in `src/tools/`. |
| `M_wormCP2AFMmu2dplaq_<tag>.dat` | Per-measurement winding numbers, charges and worm counters, comma-separated. |
| `run_<tag>.toml` | The exact settings, including the seed, that reproduce the run. |

`<tag>` encodes the parameters, e.g.
`BC_closed_nMeas_100000_beta_2.50_XYT_4_4_200_mu3_0.000_mu8_0.000_nSim_1`, and the analysis
tools read the lattice size and other values back out of it, so do not rename the files.
A run overwrites existing files with the same name; use `--output-dir` or a different
`sim` to keep earlier results.

`ctest` runs a fixed-seed case and compares it byte for byte with `tests/regression/golden`,
checks that a recorded run replays exactly, that several runs in one file match separate
runs, and that bad configs are rejected with a clear message.

### Analysing M_ output
```zsh
cd analysis_with_python && uv sync
uv run python analyse_m.py ../data/exact     # M_*.dat files or directories
```
With `q3`, `q8` the temporal winding charges (doubled units), `Lx` the periodic length, `N = Lx*Ly`:
`M3 = <q3>/(2 Lx)`, `chi3 = beta <q3^2>/(4 N)`, `M8 = sqrt(3) <q8>/(2 Lx)`, `chi8 = 3 beta <q8^2>/(4 N)`.
Errors come from a 50-block jackknife (`--nblocks`), the first 1% of rows is discarded (`--skip`).
It writes `m_summary.csv`, `susceptibility.png` and `magnetisation.png` (MC and exact values from
`tests/exact/exactCP2_results.csv` on the same axes) to `--out` (default `data/exact/plots`).

`plot_m_raw.py` draws the raw columns for checking thermalisation and autocorrelation by eye: each
panel is one column, with every `--block` (default 100) measurements averaged into one point and the
whole-run mean dashed. A drift at the start suggests thermalisation, slow wander suggests autocorrelation.
```zsh
uv run python plot_m_raw.py ../data/exact    # writes raw_<lattice>_<bc>.png to data/exact/plots
```

`plot_m_acf.py` draws the autocorrelation function of each column and estimates the integrated
autocorrelation time `tau_int = 1/2 + sum rho(t)` with Sokal's automatic window (dashed line);
it writes `acf_<lattice>_<bc>.png` and `tau_int.csv`. Blocks used for errors should be much longer than `2 tau_int`.
```zsh
uv run python plot_m_acf.py ../data/exact [--maxlag 50]
```

`thermalisation.py` helps choose how many initial measurements to skip: `skipscan_<lattice>_<bc>.png` shows
chi3 and M3 against the number skipped (they should be flat once past the transient, exact value dashed), and
`mser.csv` gives the MSER-5 truncation point per column. MSER is noisy for rare-spike columns such as `w3x`,
so judge by the skip scan and the q3/q8 columns.
```zsh
uv run python thermalisation.py ../data/exact
```

### Notes
- Boost 1.74 needs `_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION` under C++17 with libc++
  (set in `src/CMakeLists.txt`); it can be dropped after upgrading Boost.
- If the compiler cannot find `<iostream>` on macOS, a stale
  `/Library/Developer/CommandLineTools/usr/include/c++/v1` is shadowing the SDK headers.
  Reinstall the Command Line Tools, or pass `-nostdinc++ -isystem $(xcrun --show-sdk-path)/usr/include/c++/v1`.
