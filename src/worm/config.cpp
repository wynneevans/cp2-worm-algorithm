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
  if (r.sim < 0) throw InputError(where + ": sim must not be negative");
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

  rejectUnknownKeys(root, {"defaults", "run"}, path);
  const toml::table* defaults = nullptr;
  if (const toml::node* d = root.get("defaults")) {
    defaults = d->as_table();
    if (!defaults) throw InputError(path + ": [defaults] must be a table");
    rejectUnknownKeys(*defaults, {"bc", "measurements", "probabilities", "output", "mu3", "mu8", "sim", "seed"},
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
    rejectUnknownKeys(*t, {"name", "bc", "measurements", "beta", "lattice", "probabilities", "output", "mu3", "mu8",
                           "sim", "seed", "enabled"},
                      where);
    if (!names.insert(r.name).second) throw InputError(where + ": name is used by more than one run");

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
    r.sim = static_cast<int>(fetchWithDefault<int64_t>(*t, defaults, "sim", where).value_or(1));
    r.enabled = fetch<bool>(*t, "enabled", where).value_or(true);
    if (auto s = fetchWithDefault<int64_t>(*t, defaults, "seed", where)) {
      if (*s < 1 || *s > static_cast<int64_t>(kMaxSeed)) throw InputError(where + ": seed must be between 1 and 2147483647");
      r.seed = static_cast<unsigned>(*s);
    }

    const toml::array* lat = (*t)["lattice"].as_array();
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

    r.probPath = resolvePath(base, need(fetchWithDefault<std::string>(*t, defaults, "probabilities", where), "probabilities"));
    r.outDir = resolvePath(base, need(fetchWithDefault<std::string>(*t, defaults, "output", where), "output"));

    validateRun(r);
    out.push_back(r);
  }
  return out;
}

std::string datasetTag(const RunConfig& r) {
  char buf[256];
  std::snprintf(buf, sizeof buf, "BC_%s_nMeas_%d_beta_%.2f_XYT_%d_%d_%d_mu3_%.3f_mu8_%.3f_nSim_%d", r.bc.c_str(),
                r.measurements, r.beta, r.Lx, r.Ly, r.Lt, r.mu3, r.mu8, r.sim);
  return buf;
}

std::string toToml(const RunConfig& r) {
  std::ostringstream o;
  o << "# Reproduce with: wormCP2Mu2dAFMplaq <this file>\n"
    << "[[run]]\n"
    << "name = " << toml::value<std::string>(r.name) << "\n"
    << "bc = " << toml::value<std::string>(r.bc) << "\n"
    << "measurements = " << r.measurements << "\n"
    << "beta = " << formatDouble(r.beta) << "\n"
    << "lattice = [" << r.Lx << ", " << r.Ly << ", " << r.Lt << "]\n"
    << "mu3 = " << formatDouble(r.mu3) << "\n"
    << "mu8 = " << formatDouble(r.mu8) << "\n"
    << "sim = " << r.sim << "\n";
  if (r.seed) o << "seed = " << *r.seed << "\n";
  o << "probabilities = " << toml::value<std::string>(fs::absolute(r.probPath).string()) << "\n"
    << "output = " << toml::value<std::string>(fs::absolute(r.outDir).string()) << "\n";
  return o.str();
}
