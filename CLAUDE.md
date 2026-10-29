# cp2-worm-algorithm-cpp

A worm-algorithm Monte Carlo simulator for the CP(2) spin system in C++20 (`wormCP2Mu2dAFMplaq`).
The Python package `cp2analysis` (CLI `cp2-analysis`, in `analysis/`) analyses its output.
Legacy C++ analysis tools live in `src/tools/`.

## Layout
- `src/worm/` - the simulator: `main.cpp` (options, planning the runs), `runner.*` (running, `--jobs`, publishing the files), `lattice.*` (the algorithm, legacy code; `simulate()` is its only entry point), `config.*` (TOML run files), `probabilities.*` (probability table)
- `src/tools/` - legacy analysis tools; they parse the old long file names and cannot read the new run folders
- `include/cp2worm/` - shared headers
- `third_party/` - vendored ranlxd and toml++. **Do not edit.**
- `configs/` - run configs (`*.toml`) and `probabilities.txt`, the single master probability table
- `data/` - simulation output (gitignored), one folder per parameter set: `data/<bc>_<nmeas>_<beta>_<Lx>_<Ly>_<Lt>_mu3_<mu3>_mu8_<mu8>/`, each run adding `M_<UTC finish time>.dat`, `C_<time>.dat`, `run_<time>.toml` (and `log_<time>.txt` with `--jobs` > 1). Holds the user's results (GBs). **Never delete or overwrite anything here.**
- `tests/regression/` - CTest scripts, golden output, bad-config cases
- `tests/bench/` - `bench.toml`, fixed-seed runs timed by `scripts/bench.py` (not in ctest); `results/` holds the saved timings with full parameters (committed; add one per refactoring step)
- `tests/exact/` - exact-results CSV, and the `exact_*` ctest physics checks (`exact_runs.toml`, run through `cp2-analysis check-exact`, so they need uv)
- `analysis/` - Python package; tests in `analysis/tests/`
- `reference/` - frozen original simulator and the old analysis notebook. **Do not edit.**
- `scripts/` - `convert_probs.py`, `merge_probs.py`, `convert_mathematica_probs.py` (notebook output to table rows), `bench.py` (timing)
- `build/` - build folders. **Never commit, never edit by hand.**

## Commands

### C++ (run from the repo root)
- Configure: `cmake --preset=release-local` (machine-specific, from the gitignored `CMakeUserPresets.json`; presets are version 3, so CMake 3.21+; this machine has CMake 4 from Homebrew; the presets are configure-only)
- Build: `cmake --build build/release`
- Test: `(cd build/release && ctest --output-on-failure -j4)` (the `exact` tests take ~30 s; `-LE exact` skips them)
- Dry run a config: `build/release/src/wormCP2Mu2dAFMplaq configs/<file>.toml --dry-run`
- Benchmark: `scripts/bench.py --save <label>` (about 6 min; `--only long` runs just the 60 s case), then after the change `--compare tests/bench/results/<file>.json --save <label>`; `--list` shows the saved results
- macOS: if `<iostream>` is not found, a stale `/Library/Developer/CommandLineTools/usr/include/c++/v1` is hiding the SDK headers; reinstall the Command Line Tools (fixed on this machine on 2026-10-28, so the local preset no longer needs `-nostdinc++`).

### Python (run from the repo root)
- Setup: `uv sync --project analysis`
- Test: `uv run --project analysis pytest -q analysis`
- Lint and format: `uv run --project analysis ruff check analysis` and `... ruff format analysis`
- Run: `uv run --project analysis cp2-analysis {summary,raw,acf,thermalisation,check-exact} [inputs]`

Before saying a task is done, build and run the tests for whichever side you changed (C++: build + ctest; Python: pytest + ruff).

## C++ conventions
- Standard: C++20. No external C++ dependencies beyond `third_party/`; output lines are written with `std::format` (the `M_`/`C_` format strings must keep the files byte-identical; the golden tests check this).
- Bad user input (config, probability table, options) throws `InputError`, which `main` reports as `error: ...` with exit status 1. Validate everything before the first run starts.
- A state the algorithm should never reach (an unphysical branch, a bug) calls `Lattice::impossible("<function>: <what>")`, which throws `AlgorithmError`: the run stops with exit status 1 and its files stay `*_running_<pid>`, so it is never published. Don't print and carry on.
- New config keys need validation and a `bad_config_*` test in `tests/CMakeLists.txt`.
- The golden tests run each case in `tests/regression/golden/golden.toml` (closed 4x4; y-open with even and odd Ly; y-open at mu3 = 0.2 and at mu8 = 0.2) and compare the M_ and C_ files with `golden/{M,C}_<name>.dat` byte for byte. Two longer cases (`cover_mu3`, `cover_mu8`, 16x16x400) reach rare worm states the small ones miss (found with a coverage build); their output is too big to store, so they are checked against `golden/{M,C}_<name>.sha256`. They are the safety net for refactoring `lattice.cpp`. Keep them identical unless the physics is meant to change; then regenerate as described in `tests/regression/run_golden.cmake` and say so.
- `lattice.cpp` is legacy (globals, `using namespace`, duplicated pattern tables). Don't refactor it unasked.
- Changes must not slow the simulator down: compare with the latest result in `tests/bench/results/` (or save a new baseline first if the machine or build changed) using `scripts/bench.py --compare`, on an idle machine, and drop a change that is slower by more than the spread. Keep the update loop in one translation unit, with no new run-time branches, virtual calls or `std::function` in it.
- The random number generator (ranlxd) has global state: one chain per process. `--jobs N` forks one child process per run (`runParallel` in `runner.cpp`).

## Python conventions
- Python 3.11+, type hints on public functions.
- polars for data (not pandas), numpy for numerics.
- Plots with seaborn; matplotlib calls only inside `plotting.py`. Save to files, never `plt.show()`.
- Each command in `commands/` splits into compute (no I/O), plot and `run(args)`, so tests can call the logic without the CLI.
- Paths come from arguments; defaults are relative to the repo root.

## C++ / Python boundary
- `M_<time>.dat` columns are `w3x, w8x, w3y, w8y, q3, q8, l, lp, lm, tail` (written by `Lattice::measure` in `lattice.cpp`), mirrored by `COLUMNS` in `analysis/src/cp2analysis/runs.py`.
- Run parameters come from `run_<time>.toml` (written by `toToml` in `config.cpp`, read by `read_run` in `runs.py`).
- **If you change an output format, update both sides and both test suites in the same change.**

## Don't
- Edit `third_party/`, `reference/` or anything under `build/`.
- Delete or overwrite anything in `data/`. When testing, run into a scratch folder with `--data-dir`.
- Commit. The user commits (`git commit` is denied in the permissions).
- Add dependencies (C++ or Python) without asking.
- Change compiler flags or presets unless asked.
