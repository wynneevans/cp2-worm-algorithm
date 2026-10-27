#pragma once
#include <optional>
#include <string>
#include <vector>

// One simulation, fully resolved (defaults merged in, paths made concrete).
struct RunConfig {
  std::string name;
  std::string bc;  // "closed" or "y-open"
  int measurements = 0;
  double beta = 0;
  int Lx = 0, Ly = 0, Lt = 0;
  double mu3 = 0, mu8 = 0;
  std::string probPath;
  std::string outDir;
  int sim = 1;  // tag in the output file names
  std::optional<unsigned> seed;
  bool enabled = true;

  double dt() const { return 4.0 * beta / Lt; }
};

constexpr unsigned kMaxSeed = 2147483647u;  // RANLUX seeds are 1 .. 2^31-1

// Parse a TOML run file. Relative paths are taken relative to the file.
// Throws InputError with a message naming the file and run on any problem.
std::vector<RunConfig> loadConfig(const std::string& path);

// Checks a single run (also used after command-line overrides).
void validateRun(const RunConfig& run);

// "BC_closed_nMeas_1000_beta_2.50_XYT_4_4_200_mu3_0.000_mu8_0.000_nSim_1",
// the part of the output file names that identifies the run.
std::string datasetTag(const RunConfig& run);

// Settings that reproduce this run exactly (includes the seed actually used).
std::string toToml(const RunConfig& run);
