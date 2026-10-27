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
third_party/ranlxd Martin Luescher's RANLUX random number generator (GPL)
tests/             CTest regression tests and golden output
scripts/           run_worm.zsh, the batch of production runs
analysis/          Jupyter notebook and fitting scripts
reference/         frozen copy of the original /dev/urandom-seeded simulator (not built)
data/              simulation output (gitignored)
```

## Building and running
Requires CMake 3.19+, a C++17 compiler and Boost 1.74+ (headers only).

```zsh
cmake -S . -B build/release -DBoost_ROOT=/path/to/boost
cmake --build build/release
(cd build/release && ctest --output-on-failure)
```

Machine-specific settings (Boost location, compiler flags) belong in a
`CMakeUserPresets.json`, which is gitignored. Copy the `local-base` pattern from
`CMakePresets.json`'s `release` preset and then run `cmake --preset=release-local`.

The simulator takes the RNG seed as its 11th argument:
```zsh
build/release/src/wormCP2Mu2dAFMplaq data closed 100000 2.50 4 4 200 probabilities.txt 1 1 <seed>
```
`scripts/run_worm.zsh` passes a fresh random seed per run. `ctest` runs a fixed-seed
case (seed 42) and compares it byte for byte with `tests/regression/golden`, so any
change that alters the random stream or the physics fails the test.

### Notes
- Boost 1.74 needs `_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION` under C++17 with libc++
  (set in `src/CMakeLists.txt`); it can be dropped after upgrading Boost.
- If the compiler cannot find `<iostream>` on macOS, a stale
  `/Library/Developer/CommandLineTools/usr/include/c++/v1` is shadowing the SDK headers.
  Reinstall the Command Line Tools, or pass `-nostdinc++ -isystem $(xcrun --show-sdk-path)/usr/include/c++/v1`.
