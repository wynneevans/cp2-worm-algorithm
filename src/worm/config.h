#pragma once
#include <optional>
#include <string>
#include <vector>

#include "probabilities.h"

// One simulation, fully resolved (defaults merged in, paths made concrete).
struct RunConfig {
  std::string name;
  std::string bc;        // "closed" or "y-open"
  int measurements = 0;
  double beta = 0;
  int Lx = 0, Ly = 0, Lt = 0;
  double mu3 = 0, mu8 = 0;
  std::optional<std::string> probPath;              // master table; unset when the probabilities are inline
  std::optional<ProbabilityTable::Row> inlineRow;  // probabilities written in the config (a run record)
  std::string dataDir;                              // root of the output tree
  std::optional<unsigned> seed;
  bool enabled = true;

  double dt() const { return 4.0 * beta / Lt; }
  // Folder name for this parameter set, e.g. "closed_1000000_0.05_2_2_4_mu3_0.000_mu8_0.000".
  // Every run with these parameters adds its files to <dataDir>/<folderName()>.
  std::string folderName() const;
  std::string runDir() const;
};

// Facts about a finished (or started) run, written under [provenance] and ignored when replaying.
struct Provenance {
  std::string started;                  // UTC, ISO 8601
  std::optional<std::string> finished;  // absent until the run completes
  std::optional<double> wallSeconds;
  std::string git;
};

constexpr unsigned kMaxSeed = 2147483647u;  // RANLUX seeds are 1 .. 2^31-1

// Parse a TOML run file. Relative paths are taken relative to the file. When no
// probabilities are given, the nearest probabilities.txt in the file's folder
// or a parent folder is used. Throws InputError with a message naming the file
// and run on any problem.
std::vector<RunConfig> loadConfig(const std::string& path);

// Checks a single run (also used after command-line overrides).
void validateRun(const RunConfig& run);

// The row to simulate with: the inline one, or the table's row for this run.
// `table` may be null only if the run has inline probabilities.
const ProbabilityTable::Row& probabilityRow(const RunConfig& run, const ProbabilityTable* table);

// Text of the run record (run_<datetime>.toml): settings that reproduce this run
// exactly (the seed actually used, the probabilities used) plus provenance.
// data_dir is written relative to the run folder so the data tree can be moved.
std::string toToml(const RunConfig& run, const ProbabilityTable::Row& row, const Provenance& provenance);

// UTC time in ISO 8601, e.g. "2026-10-28T10:15:00Z".
std::string utcNow();

// UTC time for file names, e.g. "20261028-101500".
std::string utcStamp();
