#pragma once
#include <array>
#include <stdexcept>
#include <string>

#include "config.h"
#include "probabilities.h"

// Thrown when the algorithm reaches a state a correct simulation never reaches (a bug, not bad
// input). The run stops and its files keep their running_<pid> names, so it is not published.
struct AlgorithmError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// The worm algorithm (legacy code in lattice.cpp). Simulates `run` with the update weights `p` and
// writes M_<tag>.dat and C_<tag>.dat in `dir`. Seeds ranlxd from *run.seed, which must be set;
// ranlxd has global state, so a process runs one simulation at a time.
void simulate(const RunConfig& run, const std::array<double, kNumWeights>& p, const std::string& dir,
              const std::string& tag);
