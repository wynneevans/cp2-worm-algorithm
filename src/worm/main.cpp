#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "config.h"
#include "input_error.h"
#include "lattice.h"
#include "probabilities.h"
#include "runner.h"

using namespace std;

namespace {

struct Options {
  string configPath;
  std::optional<string> only;
  std::optional<string> dataDir;
  std::optional<unsigned> seed;
  unsigned jobs = 1;
  bool dryRun = false;
  bool help = false;
};

void printUsage(ostream& os) {
  os << "Usage: wormCP2Mu2dAFMplaq <run.toml> [options]\n"
        "\n"
        "Runs every enabled [[run]] in the TOML file, in order. Each run adds M_<time>.dat,\n"
        "C_<time>.dat and run_<time>.toml (time = UTC finish time) to the folder for its\n"
        "parameters, <data dir>/<bc>_<measurements>_<beta>_<Lx>_<Ly>_<Lt>_mu3_<mu3>_mu8_<mu8>/.\n"
        "\n"
        "Options:\n"
        "  --only NAME         run only the run called NAME (even if it is disabled); for a run\n"
        "                      with repeats, NAME selects all its copies and NAME_repK one copy\n"
        "  --data-dir DIR      put the run folders under DIR instead of the file's data_dir\n"
        "  --seed N            use seed N (1 to 2147483647) for every run\n"
        "  --jobs N            run up to N runs at once, each in its own process (default 1);\n"
        "                      each run's output then goes to log_<time>.txt in its folder\n"
        "  --dry-run           print the resolved runs and exit\n"
        "  -h, --help          show this help\n"
        "\n"
        "Without a seed (here or in the file) each run draws a random one, prints it and\n"
        "records it in run_<time>.toml, which reproduces the run.\n";
}

// The value of a numeric option as an integer from 1 to max, or nothing if it is not one.
std::optional<unsigned> parseCount(const string& v, unsigned long long max) {
  char* end = nullptr;
  errno = 0;
  const unsigned long long n = strtoull(v.c_str(), &end, 10);
  if (v.empty() || *end != '\0' || errno != 0 || n < 1 || n > max) return std::nullopt;
  return static_cast<unsigned>(n);
}

Options parseArgs(int argc, char** argv) {
  Options o;
  auto value = [&](int& i, const string& flag) -> string {
    if (i + 1 >= argc) throw InputError("option " + flag + " needs a value");
    return argv[++i];
  };
  for (int i = 1; i < argc; i++) {
    const string a = argv[i];
    if (a == "-h" || a == "--help") {
      o.help = true;
    } else if (a == "--dry-run") {
      o.dryRun = true;
    } else if (a == "--only") {
      o.only = value(i, a);
    } else if (a == "--data-dir") {
      o.dataDir = value(i, a);
    } else if (a == "--seed") {
      const string v = value(i, a);
      o.seed = parseCount(v, kMaxSeed);
      if (!o.seed) throw InputError("--seed must be an integer between 1 and 2147483647, got '" + v + "'");
    } else if (a == "--jobs") {
      const string v = value(i, a);
      const auto jobs = parseCount(v, 1024);
      if (!jobs) throw InputError("--jobs must be an integer of at least 1, got '" + v + "'");
      o.jobs = *jobs;
    } else if (a.size() > 1 && a[0] == '-') {
      throw InputError("unknown option '" + a + "'");
    } else if (o.configPath.empty()) {
      o.configPath = a;
    } else {
      throw InputError("unexpected extra argument '" + a + "'");
    }
  }
  if (!o.help && o.configPath.empty()) throw InputError("no config file given");
  return o;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc == 1) {
    printUsage(cout);
    return 0;
  }
  try {
    const Options opt = parseArgs(argc, argv);
    if (opt.help) {
      printUsage(cout);
      return 0;
    }

    vector<RunConfig> runs = loadConfig(opt.configPath);
    if (opt.only) {
      vector<RunConfig> chosen;
      for (RunConfig r : runs) {
        if (r.name != *opt.only && r.group != *opt.only) continue;
        r.enabled = true;
        chosen.push_back(r);
      }
      if (chosen.empty()) throw InputError("no run named '" + *opt.only + "' in " + opt.configPath);
      runs = chosen;
    }

    // Check every selected run before starting any of them, so a mistake in a
    // later run is reported now, not after the earlier runs have finished.
    map<string, ProbabilityTable> tables;
    Plan plan;
    for (RunConfig r : runs) {
      if (!r.enabled) {
        cout << "skipping disabled run '" << r.name << "'" << endl;
        continue;
      }
      if (opt.dataDir) r.dataDir = filesystem::path(*opt.dataDir).lexically_normal().string();
      if (opt.seed) r.seed = opt.seed;
      validateRun(r);
      const ProbabilityTable* table = nullptr;
      if (!r.inlineRow) {
        auto t = tables.find(*r.probPath);
        if (t == tables.end()) t = tables.emplace(*r.probPath, loadProbabilities(*r.probPath)).first;
        table = &t->second;
      }
      plan.emplace_back(r, probabilityRow(r, table));
    }
    if (plan.empty()) throw InputError("no enabled runs in " + opt.configPath);

    // Two runs with the same parameters and the same fixed seed would give identical output.
    map<pair<string, unsigned>, string> seeded;
    for (const auto& [r, row] : plan) {
      if (!r.seed) continue;
      const auto [it, fresh] = seeded.emplace(make_pair(r.folderName(), *r.seed), r.name);
      if (!fresh) {
        throw InputError("runs '" + it->second + "' and '" + r.name + "' have the same parameters and seed " +
                         to_string(*r.seed) + ", so they would produce identical output");
      }
    }

    if (opt.jobs > 1 && !opt.dryRun) return runParallel(plan, opt.jobs);
    for (const auto& [run, row] : plan) {
      describe(run, row);
      if (!opt.dryRun) runOne(run, row);
    }
    return 0;
  } catch (const InputError& e) {
    cerr << "error: " << e.what() << endl;
    return 1;
  } catch (const AlgorithmError& e) {
    cerr << "error: " << e.what() << endl;
    return 1;
  }
}
