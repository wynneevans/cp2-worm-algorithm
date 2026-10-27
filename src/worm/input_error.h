#pragma once
#include <stdexcept>

// Thrown for anything wrong with user-supplied input: the run config, the
// probability table or the command line. main() reports it and exits non-zero.
struct InputError : std::runtime_error {
  using std::runtime_error::runtime_error;
};
