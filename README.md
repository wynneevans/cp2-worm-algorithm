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
tests/exact/       exact results (CSV) for Monte Carlo comparison; the runs are in configs/exact.toml
configs/           run configs you edit (*.toml) and probabilities.txt, the one master probability table
scripts/           convert_probs.py, converts old probability files
analysis/  Python (uv, polars, seaborn) analysis of the output files
reference/         archive: original /dev/urandom-seeded simulator (not built) and the old analysis/ notebook
data/              simulation output, one folder per parameter set (gitignored)
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
build/release/src/wormCP2Mu2dAFMplaq configs/production.toml --only closed_b2.5_4x4x200 --data-dir data/test
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
The simulator reads a TOML run file from `configs/`; `configs/production.toml` is the batch that
used to be `run_worm.zsh`:

```zsh
build/release/src/wormCP2Mu2dAFMplaq configs/production.toml            # all enabled runs
build/release/src/wormCP2Mu2dAFMplaq configs/production.toml --only closed_b5_8x8x400
build/release/src/wormCP2Mu2dAFMplaq configs/production.toml --dry-run   # show what would run
```

```toml
[defaults]
bc            = "closed"            # or "y-open"
measurements  = 100000              # at least 20
# beta, lattice                     # may also go here, for scans that share them
mu3 = 0.0
mu8 = 0.0
# data_dir = "../data"              # output root, relative to this file (the default)
# probabilities = "probabilities.txt"   # optional, see below

[[run]]
name    = "b2.5_4x4x200"            # for --only and labels; unique within the file
beta    = 2.50                      # at most 2 decimal places (mu3, mu8: at most 3)
lattice = [4, 4, 200]               # Lx, Ly, Lt (Lt must be a multiple of 4)
# seed = 42                         # optional
# repeats = 4                       # optional: 4 independent copies, named b2.5_4x4x200_rep1 ... _rep4
# enabled = false                   # skip this run unless chosen with --only
```
Results go in one folder per parameter set, `<data_dir>/<bc>_<measurements>_<beta>_<Lx>_<Ly>_<Lt>_mu3_<mu3>_mu8_<mu8>/`,
e.g. `data/closed_100000_2.50_4_4_200_mu3_0.000_mu8_0.000/`. The folder is created by the first run
with those parameters; later runs with the same parameters add their files to it. Because the name
rounds beta to 2 and mu3, mu8 to 3 decimal places, a config with more precision is rejected.

### Command-line options
```
wormCP2Mu2dAFMplaq <run.toml> [options]
```
| Option | Effect |
| --- | --- |
| `--only NAME` | Run only the `[[run]]` called `NAME`, even if it has `enabled = false`. For a run with `repeats`, `NAME` selects all its copies and `NAME_repK` one copy. |
| `--data-dir DIR` | Put the parameter folders under `DIR` instead of the file's `data_dir` (e.g. a scratch folder for tests). |
| `--seed N` | Use seed `N` (1 to 2147483647) for every run. |
| `--jobs N` | Run up to `N` runs at once, each in its own process (default 1, one after another). |
| `--dry-run` | Print the resolved runs, their folders and the probability row each will use, then exit. |
| `-h`, `--help` | Show the usage text. |

The program exits with status 1 and an `error:` message on bad input (config, probability
table, existing run folder or options), before any simulation starts.

**Using several cores.** One run uses one core. `--jobs 4` starts the first 4 runs of the file at once and
starts the next as soon as one finishes, until all are done; the output is the same as running them one
after another. Each run's terminal output (seed, progress bar, duration) then goes to `log_<time>.txt` in
its folder (`log_running_<pid>.txt` while it runs, so `tail -f` it to follow progress), and the terminal
shows one line per run started, finished or failed. The exit status is 1 if any run failed. On this
machine (4 performance cores) 4 is a good choice.

**Repeats.** `repeats = N` in a `[[run]]` makes N independent runs of that parameter set, named
`<name>_rep1` ... `<name>_repN`, each with its own random seed; they share one parameter folder, and
`cp2-analysis` pools them. Two runs with the same parameters and the same fixed seed (for example
`repeats` together with a `seed`, or `--seed`) would give identical output, so that is rejected.

Any `[defaults]` key can be overridden inside a `[[run]]` (`repeats` and `enabled` belong to a run only). Unknown keys are errors, so a
typo cannot silently fall back to a default. Everything is checked before the first run
starts.

**Probabilities.** `configs/probabilities.txt` is the single master table, with a header line and one
row per `(dt, mu3, mu8)`, where `dt = 4*beta/Lt`:
```
# comments start with '#'
dt mu3 mu8 p01 p02 ... p48
0.05 0 0 1. 0 0.0253044 ...
```
A run finds its row from its own `4*beta/Lt`, `mu3` and `mu8`, and stops with the list of available
rows if there is none. Without a `probabilities` key the nearest `probabilities.txt` in the config's
folder or a parent folder is used. To simulate a new `dt`, add a row. Older headerless files convert
with `scripts/convert_probs.py old.txt new.txt`; `scripts/merge_probs.py` merges tables.
A run may also give its probabilities inline, as a `[run.probabilities]` table with `dt`, `mu3`, `mu8`
and `p01` ... `p48` (this is what a run record contains).

### Output
Each finished run adds three files to its parameter folder, named by the UTC time it finished:

| File | Contents |
| --- | --- |
| `M_<time>.dat` | Per-measurement winding numbers, charges and worm counters, comma-separated. |
| `C_<time>.dat` | Correlator measurements, comma-separated. |
| `log_<time>.txt` | Only with `--jobs` above 1: the run's terminal output. |
| `run_<time>.toml` | A record of what ran, not a config you edit: the settings, the seed actually used, the probabilities used (inline) and a `[provenance]` table (start and finish time, wall seconds, git description from when CMake last configured). |

`<time>` looks like `20261028-114625`. If another run with the same parameters finished in the same
second, `-2`, `-3`, ... is added; files are never overwritten. While a run is going its files are named
`M_running_<process id>.dat` etc., so a crashed or killed run leaves `running` files behind rather than
something that looks finished.

**Replay.** A run record is itself a valid run file and reproduces the run byte for byte:
`wormCP2Mu2dAFMplaq data/<folder>/run_<time>.toml`. Without `--data-dir` the replay lands in the same
folder as a new set of files. A run with no `seed` draws a random one, prints it and records it.

The old flat file names (`M_wormCP2AFMmu2dplaq_BC_..._nSim_1.dat`) are no longer written; the 2020–21
results with those names are in `data/old_macbook_run/`. Only the C++ tools in `src/tools` read them
(they parse the names, so they cannot read the new folders); `cp2-analysis` reads the new layout only.

`ctest` runs a fixed-seed case and compares it byte for byte with `tests/regression/golden`,
checks that a recorded run replays exactly and that repeated runs with the same parameters share a
folder without overwriting each other, that
several runs in one file match separate runs, that the master table is found in a parent folder, and
that bad configs are rejected with a clear message. The tests labelled `exact` (about 30 s in all;
`ctest -j4` runs them in parallel, `ctest -LE exact` skips them) compare tiny lattices with the exact
results in `tests/exact/` using `cp2-analysis check-exact`, so they need uv (they show as disabled
without it); see `tests/exact/README.md`.

### Analysing M_ output
The `cp2-analysis` command (package in `analysis/`, managed with uv) analyses the
`M_<time>.dat` files. Run it from the repo root. Inputs default to `data/`; choose runs by folder or by
parameter (`--beta 56.3 --mu3 0.01:0.03`). Repeat runs of a parameter set are pooled unless `--per-run`:
```zsh
uv sync --project analysis
uv run --project analysis cp2-analysis summary          # magnetisation and susceptibility vs exact
uv run --project analysis cp2-analysis summary --beta 56.3   # a scan, plotted against mu3
uv run --project analysis cp2-analysis raw              # block-averaged raw columns
uv run --project analysis cp2-analysis acf              # autocorrelation function, tau_int
uv run --project analysis cp2-analysis thermalisation   # skip scan and MSER
uv run --project analysis cp2-analysis check-exact --beta 0.05 --typed   # pass/fail vs exact results
```
See [analysis/README.md](analysis/README.md) for the estimators, options and outputs.

### Notes
- Boost 1.74 needs `_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION` under C++17 with libc++
  (set in `src/CMakeLists.txt`); it can be dropped after upgrading Boost.
- If the compiler cannot find `<iostream>` on macOS, a stale
  `/Library/Developer/CommandLineTools/usr/include/c++/v1` is shadowing the SDK headers.
  Reinstall the Command Line Tools, or pass `-nostdinc++ -isystem $(xcrun --show-sdk-path)/usr/include/c++/v1`.
