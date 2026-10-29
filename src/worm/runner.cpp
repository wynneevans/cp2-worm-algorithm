#include "runner.h"

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>

#include "input_error.h"
#include "lattice.h"

using namespace std;

string runningTag() { return "running_" + to_string(getpid()); }

namespace {

unsigned randomSeed() {
  random_device rd;
  return rd() % kMaxSeed + 1;
}

// Gives the finished run's files their final names, M_<stem>.dat, C_<stem>.dat,
// run_<stem>.toml and, for a run started with --jobs, log_<stem>.txt, where stem is the UTC finish time, with -2, -3, ... added if another
// run already used that second. Uses link(), which never replaces an existing file, so
// runs finishing at the same moment in different processes cannot overwrite each other.
string publish(const filesystem::path& dir, const string& tmp) {
  const string base = utcStamp();
  auto name = [&](const string& kind, const string& stem) {
    const string ext = kind == "run" ? ".toml" : kind == "log" ? ".txt" : ".dat";
    return (dir / (kind + "_" + stem + ext)).string();
  };
  vector<string> kinds = {"C", "run"};
  if (filesystem::exists(name("log", tmp))) kinds.push_back("log");
  for (int n = 1;; n++) {
    const string stem = n == 1 ? base : base + "-" + to_string(n);
    if (link(name("M", tmp).c_str(), name("M", stem).c_str()) != 0) {
      if (errno == EEXIST) continue;
      throw InputError("cannot create '" + name("M", stem) + "': " + strerror(errno));
    }
    for (const string& kind : kinds) {
      if (link(name(kind, tmp).c_str(), name(kind, stem).c_str()) != 0) {
        throw InputError("cannot create '" + name(kind, stem) + "': " + strerror(errno));
      }
    }
    unlink(name("M", tmp).c_str());
    for (const string& kind : kinds) unlink(name(kind, tmp).c_str());
    return stem;
  }
}

}  // namespace

void describe(const RunConfig& r, const ProbabilityTable::Row& row) {
  cout << "run '" << r.name << "': bc=" << r.bc << " beta=" << r.beta << " lattice=" << r.Lx << "x" << r.Ly << "x"
       << r.Lt << " dt=" << r.dt() << " mu3=" << r.mu3 << " mu8=" << r.mu8 << " measurements=" << r.measurements
       << " seed=";
  if (r.seed) cout << *r.seed; else cout << "(random)";
  cout << "\n  probabilities: ";
  if (r.inlineRow) cout << "inline"; else cout << *r.probPath << " line " << row.line;
  cout << "\n  folder: " << r.runDir() << endl;
}

void runOne(RunConfig run, const ProbabilityTable::Row& row) {
  if (!run.seed) run.seed = randomSeed();
  const filesystem::path dir = run.runDir();
  filesystem::create_directories(dir);

  const filesystem::path record = dir / ("run_" + runningTag() + ".toml");
  Provenance prov{utcNow(), std::nullopt, std::nullopt, CP2WORM_GIT_DESCRIBE};
  {
    ofstream out(record);
    out << toToml(run, row, prov);
    if (!out) throw InputError("cannot write '" + record.string() + "'");
  }
  cout << "seed " << *run.seed << ", writing to " << dir.string() << endl;

  const auto start = chrono::steady_clock::now();
  simulate(run, row.p, dir.string(), runningTag());
  const double wall = chrono::duration<double>(chrono::steady_clock::now() - start).count();
  cout << "Duration: " << wall << endl;

  prov.finished = utcNow();
  prov.wallSeconds = wall;
  {
    ofstream out(record);
    out << toToml(run, row, prov);
    if (!out) throw InputError("cannot write '" + record.string() + "'");
  }
  const string stem = publish(dir, runningTag());
  cout << "results: " << (dir / ("M_" + stem + ".dat")).string() << ", C_" << stem << ".dat, run_" << stem << ".toml"
       << endl;
}

int runParallel(const Plan& plan, unsigned jobs) {
  map<pid_t, string> running;
  size_t next = 0;
  int failed = 0;
  bool forkFailed = false;
  while ((next < plan.size() && !forkFailed) || !running.empty()) {
    if (next < plan.size() && !forkFailed && running.size() < jobs) {
      const auto& [run, row] = plan[next++];
      const filesystem::path dir = run.runDir();
      describe(run, row);
      cout.flush();
      fflush(nullptr);
      const pid_t pid = fork();
      if (pid < 0) {
        cerr << "error: cannot start run '" << run.name << "': " << strerror(errno)
             << "; waiting for the running ones" << endl;
        forkFailed = true;
        failed++;
        continue;
      }
      if (pid == 0) {
        int code = 0;
        try {
          filesystem::create_directories(dir);
          const string log = (dir / ("log_" + runningTag() + ".txt")).string();
          const int fd = open(log.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
          if (fd < 0) throw InputError("cannot create '" + log + "': " + strerror(errno));
          dup2(fd, STDOUT_FILENO);
          dup2(fd, STDERR_FILENO);
          close(fd);
          runOne(run, row);
        } catch (const std::exception& e) {
          cerr << "error: " << e.what() << endl;
          code = 1;
        }
        cout.flush();
        fflush(nullptr);
        _exit(code);
      }
      running[pid] = run.name;
      cout << "started '" << run.name << "' (pid " << pid << "), log: "
           << (dir / ("log_running_" + to_string(pid) + ".txt")).string() << endl;
      continue;
    }
    int status = 0;
    const pid_t pid = waitpid(-1, &status, 0);
    if (pid < 0) {
      if (errno == EINTR) continue;
      cerr << "error: waitpid: " << strerror(errno) << endl;
      return 1;
    }
    const auto it = running.find(pid);
    if (it == running.end()) continue;
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
      cout << "finished '" << it->second << "'" << endl;
    } else {
      failed++;
      cout << "run '" << it->second << "' failed (";
      if (WIFSIGNALED(status)) cout << "signal " << WTERMSIG(status); else cout << "exit status " << WEXITSTATUS(status);
      cout << "), see its log" << endl;
    }
    running.erase(it);
  }
  if (failed) cerr << "error: " << failed << " of " << plan.size() << " runs did not finish" << endl;
  return failed ? 1 : 0;
}
