#pragma once
#include <array>
#include <string>
#include <vector>

constexpr int kNumWeights = 48;

// Table of update-probability weights, one row per (dt, mu3, mu8) where
// dt = 4*beta/Lt. File format: '#' comments, a header line
//   dt mu3 mu8 p01 p02 ... p48
// and one whitespace-separated row of numbers per parameter set.
struct ProbabilityTable {
  struct Row {
    double dt, mu3, mu8;
    std::array<double, kNumWeights> p;
    int line;  // line number in the file, for messages
  };
  std::string path;
  std::vector<Row> rows;
};

// Equal to a relative tolerance of 1e-9 (the comparison used for table lookups).
bool sameValue(double a, double b);

ProbabilityTable loadProbabilities(const std::string& path);

// The unique row matching (dt, mu3, mu8) to a relative tolerance of 1e-9.
// Throws InputError if there is none or more than one.
const ProbabilityTable::Row& selectRow(const ProbabilityTable& table, double dt, double mu3, double mu8);
