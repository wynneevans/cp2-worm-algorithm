#include "config.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <initializer_list>
#include <set>
#include <sstream>
#include <string_view>
#include <ctime>

#include <toml.hpp>

#include "input_error.h"

namespace fs = std::filesystem;

namespace {

void rejectUnknownKeys(const toml::table& t, std::initializer_list<std::string_view> allowed,
                       const std::string& where) {
  for (auto&& [key, node] : t) {
    if (std::find(allowed.begin(), allowed.end(), key.str()) == allowed.end()) {
      throw InputError(where + ": unknown key '" + std::string(key.str()) + "'");
    }
  }
}

template <class T>
std::optional<T> fetch(const toml::table& t, std::string_view key, const std::string& where) {
  const toml::node* n = t.get(key);
  if (!n) return std::nullopt;
  if (auto v = n->value<T>()) return v;
  throw InputError(where + ": key '" + std::string(key) + "' has the wrong type or value");
}

// run's own value, else the [defaults] value.
template <class T>
std::optional<T> fetchWithDefault(const toml::table& run, const toml::table* defaults, std::string_view key,
                                  const std::string& where) {
  if (auto v = fetch<T>(run, key, where)) return v;
  if (defaults) return fetch<T>(*defaults, key, "[defaults]");
  return std::nullopt;
}

int toInt(int64_t v, const std::string& what, const std::string& where) {
  if (v < 1 || v > 1'000'000'000) throw InputError(where + ": " + what + " must be a positive integer");
  return static_cast<int>(v);
}

std::string resolvePath(const fs::path& base, const std::string& p) {
  fs::path path(p);
  if (path.is_relative()) path = base / path;
  return path.lexically_normal().string();
}

std::string formatDouble(double v) {
  char buf[40];
  auto res = std::to_chars(buf, buf + sizeof buf, v);
  std::string s(buf, res.ptr);
  if (s.find_first_of(".eEn") == std::string::npos) s += ".0";  // keep it a TOML float
  return s;
}

}  // namespace


std::string RunConfig::folderName() const {
  // + 0.0 turns -0.0 into 0.0 so it prints as "0.000".
  char buf[256];
  std::snprintf(buf, sizeof buf, "%s_%d_%.2f_%d_%d_%d_mu3_%.3f_mu8_%.3f", bc.c_str(), measurements, beta + 0.0, Lx, Ly,
                Lt, mu3 + 0.0, mu8 + 0.0);
  return buf;
}

std::string RunConfig::runDir() const { return (fs::path(dataDir) / folderName()).lexically_normal().string(); }

std::string utcStamp() {
  const std::time_t now = std::time(nullptr);
  std::tm tm{};
  gmtime_r(&now, &tm);
  char buf[32];
  std::strftime(buf, sizeof buf, "%Y%m%d-%H%M%S", &tm);
  return buf;
}

std::string utcNow() {
  const std::time_t now = std::time(nullptr);
  std::tm tm{};
  gmtime_r(&now, &tm);
  char buf[32];
  std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", &tm);
  return buf;
}

namespace {

// True if v has at most `places` decimal places (as written in the config).
bool hasAtMostDecimals(double v, int places) {
  const double scaled = v * std::pow(10.0, places);
  return std::fabs(scaled - std::round(scaled)) <= 1e-9 * std::max(1.0, std::fabs(scaled));
}

std::string weightName(int i) {
  char name[8];
  std::snprintf(name, sizeof name, "p%02d", i);
  return name;
}

// Inline probabilities: a table with dt, mu3, mu8 and p01 .. p48.
ProbabilityTable::Row parseInlineRow(const toml::table& t, const std::string& where) {
  std::set<std::string> allowed = {"dt", "mu3", "mu8"};
  for (int i = 1; i <= kNumWeights; i++) allowed.insert(weightName(i));
  for (auto&& [key, node] : t) {
    if (!allowed.count(std::string(key.str()))) {
      throw InputError(where + ": probabilities: unknown key '" + std::string(key.str()) + "'");
    }
  }
  auto get = [&](const std::string& key) {
    auto v = fetch<double>(t, key, where + ": probabilities");
    if (!v) throw InputError(where + ": probabilities: missing '" + key + "'");
    if (!std::isfinite(*v)) throw InputError(where + ": probabilities: '" + key + "' is not finite");
    return *v;
  };
  ProbabilityTable::Row row;
  row.dt = get("dt");
  row.mu3 = get("mu3");
  row.mu8 = get("mu8");
  for (int i = 0; i < kNumWeights; i++) row.p[i] = get(weightName(i + 1));
  row.line = 0;
  return row;
}

// The nearest probabilities.txt in `dir` or a parent folder.
std::optional<std::string> findMasterTable(const fs::path& dir) {
  for (fs::path d = fs::absolute(dir).lexically_normal(); !d.empty(); d = d.parent_path()) {
    if (fs::exists(d / "probabilities.txt")) return (d / "probabilities.txt").string();
    if (d == d.parent_path()) break;
  }
  return std::nullopt;
}

}  // namespace

const ProbabilityTable::Row& probabilityRow(const RunConfig& run, const ProbabilityTable* table) {
  if (run.inlineRow) {
    const auto& r = *run.inlineRow;
    if (!sameValue(r.dt, run.dt()) || !sameValue(r.mu3, run.mu3) || !sameValue(r.mu8, run.mu8)) {
      std::ostringstream msg;
      msg << "run '" << run.name << "': inline probabilities are for dt=" << r.dt << " mu3=" << r.mu3 << " mu8=" << r.mu8
          << " but the run has dt=" << run.dt() << " mu3=" << run.mu3 << " mu8=" << run.mu8;
      throw InputError(msg.str());
    }
    return r;
  }
  return selectRow(*table, run.dt(), run.mu3, run.mu8);
}

void validateRun(const RunConfig& r) {
  const std::string where = "run '" + r.name + "'";
  if (r.bc != "closed" && r.bc != "y-open") {
    throw InputError(where + ": bc must be \"closed\" or \"y-open\", got \"" + r.bc + "\"");
  }
  if (r.measurements < 1) throw InputError(where + ": measurements must be positive");
  if (!(r.beta > 0) || !std::isfinite(r.beta)) throw InputError(where + ": beta must be a positive number");
  if (r.Lx < 1 || r.Ly < 1 || r.Lt < 1) throw InputError(where + ": lattice sizes must be positive");
  if (r.Lt % 4 != 0) {
    throw InputError(where + ": Lt must be a multiple of 4 (the update patterns repeat every 4 time slices), got " +
                     std::to_string(r.Lt));
  }
  if (!std::isfinite(r.mu3) || !std::isfinite(r.mu8)) throw InputError(where + ": mu3 and mu8 must be finite");
  // The folder name rounds beta to 2 and mu to 3 decimals; more precision would mix different runs.
  if (!hasAtMostDecimals(r.beta, 2)) {
    throw InputError(where + ": beta may have at most 2 decimal places (it names the output folder), got " +
                     formatDouble(r.beta));
  }
  if (!hasAtMostDecimals(r.mu3, 3) || !hasAtMostDecimals(r.mu8, 3)) {
    throw InputError(where + ": mu3 and mu8 may have at most 3 decimal places (they name the output folder), got mu3=" +
                     formatDouble(r.mu3) + " mu8=" + formatDouble(r.mu8));
  }
  if (r.seed && (*r.seed < 1 || *r.seed > kMaxSeed)) {
    throw InputError(where + ": seed must be between 1 and 2147483647");
  }
}

std::vector<RunConfig> loadConfig(const std::string& path) {
  toml::table root;
  try {
    root = toml::parse_file(path);
  } catch (const toml::parse_error& e) {
    std::ostringstream msg;
    msg << "cannot read config: " << e.description() << " (" << e.source() << ")";
    throw InputError(msg.str());
  }
  const fs::path base = fs::path(path).parent_path();

  rejectUnknownKeys(root, {"defaults", "run", "provenance"}, path);
  const toml::table* defaults = nullptr;
  if (const toml::node* d = root.get("defaults")) {
    defaults = d->as_table();
    if (!defaults) throw InputError(path + ": [defaults] must be a table");
    rejectUnknownKeys(*defaults,
                      {"bc", "measurements", "beta", "lattice", "probabilities", "data_dir", "mu3", "mu8", "seed"},
                      "[defaults]");
  }

  const toml::array* runs = root["run"].as_array();
  if (!runs || runs->empty()) throw InputError(path + ": no [[run]] tables found");

  std::vector<RunConfig> out;
  std::set<std::string> names;
  for (size_t i = 0; i < runs->size(); i++) {
    const toml::table* t = runs->get(i)->as_table();
    if (!t) throw InputError(path + ": every 'run' entry must be a [[run]] table");

    RunConfig r;
    r.name = fetch<std::string>(*t, "name", "run #" + std::to_string(i + 1)).value_or("run" + std::to_string(i + 1));
    const std::string where = "run '" + r.name + "'";
    rejectUnknownKeys(*t, {"name", "bc", "measurements", "beta", "lattice", "probabilities", "data_dir", "mu3", "mu8",
                           "seed", "enabled", "repeats"},
                      where);
    r.group = r.name;
    const int64_t repeats = fetch<int64_t>(*t, "repeats", where).value_or(1);
    if (repeats < 1) throw InputError(where + ": repeats must be at least 1");

    auto need = [&](auto value, const char* key) {
      if (!value) throw InputError(where + ": missing required key '" + key + "'");
      return *value;
    };

    r.bc = need(fetchWithDefault<std::string>(*t, defaults, "bc", where), "bc");
    r.measurements = toInt(need(fetchWithDefault<int64_t>(*t, defaults, "measurements", where), "measurements"),
                           "measurements", where);
    r.beta = need(fetchWithDefault<double>(*t, defaults, "beta", where), "beta");
    r.mu3 = fetchWithDefault<double>(*t, defaults, "mu3", where).value_or(0.0);
    r.mu8 = fetchWithDefault<double>(*t, defaults, "mu8", where).value_or(0.0);
    r.enabled = fetch<bool>(*t, "enabled", where).value_or(true);
    if (auto s = fetchWithDefault<int64_t>(*t, defaults, "seed", where)) {
      if (*s < 1 || *s > static_cast<int64_t>(kMaxSeed)) throw InputError(where + ": seed must be between 1 and 2147483647");
      r.seed = static_cast<unsigned>(*s);
    }

    const toml::node* latNode = t->get("lattice");
    if (!latNode && defaults) latNode = defaults->get("lattice");
    if (!latNode) throw InputError(where + ": missing required key 'lattice'");
    const toml::array* lat = latNode->as_array();
    if (!lat || lat->size() != 3) throw InputError(where + ": 'lattice' must be an array [Lx, Ly, Lt]");
    int dims[3];
    for (size_t k = 0; k < 3; k++) {
      auto v = lat->get(k)->value<int64_t>();
      if (!v) throw InputError(where + ": 'lattice' entries must be integers");
      dims[k] = toInt(*v, "lattice size", where);
    }
    r.Lx = dims[0];
    r.Ly = dims[1];
    r.Lt = dims[2];

    r.dataDir = resolvePath(base, fetchWithDefault<std::string>(*t, defaults, "data_dir", where).value_or("../data"));

    const toml::node* probs = t->get("probabilities");
    if (probs && probs->is_table()) {
      r.inlineRow = parseInlineRow(*probs->as_table(), where);
    } else if (auto p = fetchWithDefault<std::string>(*t, defaults, "probabilities", where)) {
      r.probPath = resolvePath(base, *p);
    } else if (auto found = findMasterTable(base.empty() ? fs::path(".") : base)) {
      r.probPath = *found;
    } else {
      throw InputError(where + ": no probabilities: set 'probabilities' or put a probabilities.txt in the config's folder or a parent folder");
    }

    validateRun(r);
    for (int64_t k = 1; k <= repeats; k++) {
      RunConfig copy = r;
      if (repeats > 1) copy.name = r.name + "_rep" + std::to_string(k);
      if (!names.insert(copy.name).second) {
        throw InputError("run '" + copy.name + "': name is used by more than one run");
      }
      out.push_back(copy);
    }
  }
  return out;
}

std::string toToml(const RunConfig& r, const ProbabilityTable::Row& row, const Provenance& prov) {
  const fs::path runDir = fs::absolute(r.runDir()).lexically_normal();
  const fs::path dataDir = fs::absolute(r.dataDir).lexically_normal();
  std::ostringstream o;
  o << "# Record of what was run. Replay with: wormCP2Mu2dAFMplaq <this file> [--data-dir <elsewhere>]\n"
    << "[[run]]\n"
    << "name = " << toml::value<std::string>(r.name) << "\n"
    << "data_dir = " << toml::value<std::string>(dataDir.lexically_relative(runDir).generic_string()) << "\n"
    << "bc = " << toml::value<std::string>(r.bc) << "\n"
    << "measurements = " << r.measurements << "\n"
    << "beta = " << formatDouble(r.beta) << "\n"
    << "lattice = [" << r.Lx << ", " << r.Ly << ", " << r.Lt << "]\n"
    << "mu3 = " << formatDouble(r.mu3) << "\n"
    << "mu8 = " << formatDouble(r.mu8) << "\n";
  if (r.seed) o << "seed = " << *r.seed << "\n";
  o << "\n[run.probabilities]\n"
    << "dt = " << formatDouble(row.dt) << "\n"
    << "mu3 = " << formatDouble(row.mu3) << "\n"
    << "mu8 = " << formatDouble(row.mu8) << "\n";
  for (int i = 0; i < kNumWeights; i++) o << weightName(i + 1) << " = " << formatDouble(row.p[i]) << "\n";
  o << "\n[provenance]\n"
    << "started = " << toml::value<std::string>(prov.started) << "\n";
  if (prov.finished) o << "finished = " << toml::value<std::string>(*prov.finished) << "\n";
  if (prov.wallSeconds) o << "wall_seconds = " << formatDouble(*prov.wallSeconds) << "\n";
  o << "git = " << toml::value<std::string>(prov.git) << "\n";
  return o.str();
}
