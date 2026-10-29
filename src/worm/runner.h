#pragma once
#include <string>
#include <utility>
#include <vector>

#include "config.h"
#include "probabilities.h"

// Every run to do, with its probability row, resolved and checked before the first one starts.
using Plan = std::vector<std::pair<RunConfig, ProbabilityTable::Row>>;

// While a run is going its files are M_<tag>.dat, C_<tag>.dat and run_<tag>.toml with this
// tag; when it finishes they are renamed to the finish time.
std::string runningTag();

// Prints the resolved settings of a run.
void describe(const RunConfig& r, const ProbabilityTable::Row& row);

// Runs one simulation in this process (drawing a seed if it has none), writes its run record and
// gives its files their final names.
void runOne(RunConfig run, const ProbabilityTable::Row& row);

// Runs each planned run in its own process, at most `jobs` at a time, starting the next as
// soon as one finishes. Each child is a fresh single-threaded copy, so its output is the
// same as in a sequential run. A child's stdout and stderr go to log_running_<pid>.txt in
// its run folder (published as log_<time>.txt). Returns the exit status for main.
int runParallel(const Plan& plan, unsigned jobs);
