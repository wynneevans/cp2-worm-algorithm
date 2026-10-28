#include "probabilities.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "input_error.h"

namespace {

std::vector<std::string> expectedHeader() {
  std::vector<std::string> h = {"dt", "mu3", "mu8"};
  for (int i = 1; i <= kNumWeights; i++) {
    char name[8];
    std::snprintf(name, sizeof name, "p%02d", i);
    h.push_back(name);
  }
  return h;
}

std::vector<std::string> tokens(const std::string& line) {
  std::istringstream in(line);
  std::vector<std::string> out;
  for (std::string tok; in >> tok;) out.push_back(tok);
  return out;
}

}  // namespace

bool sameValue(double a, double b) {
  return std::fabs(a - b) <= 1e-9 * std::max({1.0, std::fabs(a), std::fabs(b)});
}

ProbabilityTable loadProbabilities(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw InputError("cannot open probability file '" + path + "'");

  ProbabilityTable table;
  table.path = path;
  const std::vector<std::string> header = expectedHeader();
  bool sawHeader = false;

  std::string line;
  for (int lineNo = 1; std::getline(in, line); lineNo++) {
    const std::string where = path + ":" + std::to_string(lineNo) + ": ";
    line = line.substr(0, line.find('#'));
    const std::vector<std::string> tok = tokens(line);
    if (tok.empty()) continue;

    if (!sawHeader) {
      if (tok != header) {
        throw InputError(where + "expected the header 'dt mu3 mu8 p01 ... p48' (" +
                         std::to_string(header.size()) + " columns), found " + std::to_string(tok.size()) +
                         " columns starting with '" + tok[0] + "'");
      }
      sawHeader = true;
      continue;
    }

    if (tok.size() != header.size()) {
      throw InputError(where + "expected " + std::to_string(header.size()) + " values, found " +
                       std::to_string(tok.size()));
    }
    std::vector<double> v(tok.size());
    for (size_t i = 0; i < tok.size(); i++) {
      char* end = nullptr;
      v[i] = std::strtod(tok[i].c_str(), &end);
      if (end != tok[i].c_str() + tok[i].size() || !std::isfinite(v[i])) {
        throw InputError(where + "column '" + header[i] + "': '" + tok[i] + "' is not a finite number");
      }
    }
    ProbabilityTable::Row row;
    row.dt = v[0];
    row.mu3 = v[1];
    row.mu8 = v[2];
    std::copy(v.begin() + 3, v.end(), row.p.begin());
    row.line = lineNo;
    table.rows.push_back(row);
  }

  if (!sawHeader) throw InputError(path + ": no header line found");
  if (table.rows.empty()) throw InputError(path + ": header found but no data rows");
  return table;
}

const ProbabilityTable::Row& selectRow(const ProbabilityTable& table, double dt, double mu3, double mu8) {
  const ProbabilityTable::Row* found = nullptr;
  for (const auto& r : table.rows) {
    if (sameValue(r.dt, dt) && sameValue(r.mu3, mu3) && sameValue(r.mu8, mu8)) {
      if (found) {
        throw InputError(table.path + ": lines " + std::to_string(found->line) + " and " +
                         std::to_string(r.line) + " both match dt=" + std::to_string(dt));
      }
      found = &r;
    }
  }
  if (found) return *found;

  std::ostringstream msg;
  msg << table.path << ": no row for dt=" << dt << " (4*beta/Lt), mu3=" << mu3 << ", mu8=" << mu8
      << ". Rows available (dt mu3 mu8):";
  const size_t shown = std::min<size_t>(table.rows.size(), 10);
  for (size_t i = 0; i < shown; i++) {
    msg << "\n  line " << table.rows[i].line << ": " << table.rows[i].dt << " " << table.rows[i].mu3 << " "
        << table.rows[i].mu8;
  }
  if (table.rows.size() > shown) msg << "\n  ... and " << table.rows.size() - shown << " more";
  throw InputError(msg.str());
}
