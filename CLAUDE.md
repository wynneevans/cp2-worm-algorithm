# cp2-worm-algorithm-cpp

A worm-algorithm Monte Carlo simulator for the CP(2) spin system in C++17 (`wormCP2Mu2dAFMplaq`).
The Python package `cp2analysis` (CLI `cp2-analysis`, in `analysis/`) analyses its output.
Legacy C++ analysis tools live in `src/tools/`.

## Layout
- `src/worm/` - the simulator: `main.cpp` (algorithm, legacy code), `config.*` (TOML run files), `probabilities.*` (probability table)
- `src/tools/` - legacy analysis tools; they parse the old long file names and cannot read the new run folders
- `include/cp2worm/` - shared headers
- `third_party/` - vendored ranlxd and toml++. **Do not edit.**
- `configs/` - run configs (`*.toml`) and `probabilities.txt`, the single master probability table
- `data/` - simulation output (gitignored), one folder per parameter set: `data/<bc>_<nmeas>_<beta>_<Lx>_<Ly>_<Lt>_mu3_<mu3>_mu8_<mu8>/`, each run adding `M_<UTC finish time>.dat`, `C_<time>.dat`, `run_<time>.toml` (and `log_<time>.txt` with `--jobs` > 1). Holds the user's results (GBs). **Never delete or overwrite anything here.**
- `tests/regression/` - CTest scripts, golden output, bad-config cases
- `tests/exact/` - exact-results CSV, and the `exact_*` ctest physics checks (`exact_runs.toml`, run through `cp2-analysis check-exact`, so they need uv)
- `analysis/` - Python package; tests in `analysis/tests/`
- `reference/` - frozen original simulator and the old analysis notebook. **Do not edit.**
- `scripts/` - `convert_probs.py`, `merge_probs.py`
- `build/` - build folders. **Never commit, never edit by hand.**

## Commands

### C++ (run from the repo root)
- Configure: `cmake --preset=release-local` (machine-specific, from the gitignored `CMakeUserPresets.json`; CMake 3.19, presets are configure-only)
- Build: `cmake --build build/release`
- Test: `(cd build/release && ctest --output-on-failure -j4)` (the `exact` tests take ~30 s; `-LE exact` skips them)
- Dry run a config: `build/release/src/wormCP2Mu2dAFMplaq configs/<file>.toml --dry-run`
- macOS: a stale `/Library/Developer/CommandLineTools/usr/include/c++/v1` hides the SDK headers; the local preset passes `-nostdinc++ -isystem $(xcrun --show-sdk-path)/usr/include/c++/v1`.

### Python (run from the repo root)
- Setup: `uv sync --project analysis`
- Test: `uv run --project analysis pytest -q analysis`
- Lint and format: `uv run --project analysis ruff check analysis` and `... ruff format analysis`
- Run: `uv run --project analysis cp2-analysis {summary,raw,acf,thermalisation,check-exact} [inputs]`

Before saying a task is done, build and run the tests for whichever side you changed (C++: build + ctest; Python: pytest + ruff).

## C++ conventions
- Standard: C++17 (C++20 clashes: `std::format` vs `boost::format` under `using namespace`).
- Bad user input (config, probability table, options) throws `InputError`, which `main` reports as `error: ...` with exit status 1. Validate everything before the first run starts.
- New config keys need validation and a `bad_config_*` test in `tests/CMakeLists.txt`.
- The golden test compares the M_ and C_ files with `tests/regression/golden/{M,C}.dat` byte for byte. Keep them identical unless the physics is meant to change; then regenerate as described in `tests/regression/run_golden.cmake` and say so.
- `main.cpp` is legacy (globals, `using namespace`, duplicated pattern tables). Don't refactor it unasked.
- The random number generator (ranlxd) has global state: one chain per process. `--jobs N` forks one child process per run (`runParallel` in `main.cpp`).

## Python conventions
- Python 3.11+, type hints on public functions.
- polars for data (not pandas), numpy for numerics.
- Plots with seaborn; matplotlib calls only inside `plotting.py`. Save to files, never `plt.show()`.
- Each command in `commands/` splits into compute (no I/O), plot and `run(args)`, so tests can call the logic without the CLI.
- Paths come from arguments; defaults are relative to the repo root.

## C++ / Python boundary
- `M_<time>.dat` columns are `w3x, w8x, w3y, w8y, q3, q8, l, lp, lm, tail` (written by `Lattice::measure` in `main.cpp`), mirrored by `COLUMNS` in `analysis/src/cp2analysis/runs.py`.
- Run parameters come from `run_<time>.toml` (written by `toToml` in `config.cpp`, read by `read_run` in `runs.py`).
- **If you change an output format, update both sides and both test suites in the same change.**

## Don't
- Edit `third_party/`, `reference/` or anything under `build/`.
- Delete or overwrite anything in `data/`. When testing, run into a scratch folder with `--data-dir`.
- Commit. The user commits (`git commit` is denied in the permissions).
- Add dependencies (C++ or Python) without asking.
- Change compiler flags or presets unless asked.
