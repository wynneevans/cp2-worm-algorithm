set(GOLDEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/regression/golden)

# A run with the golden seed must reproduce the committed .dat files exactly.
add_test(NAME golden_seed42
  COMMAND ${CMAKE_COMMAND}
    -DSIM=$<TARGET_FILE:wormCP2Mu2dAFMplaq>
    -DGOLDEN_DIR=${GOLDEN_DIR}
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/golden_seed42
    -DSEED=42 -DEXPECT_MATCH=ON
    -P ${CMAKE_CURRENT_SOURCE_DIR}/regression/run_golden.cmake)

# Guard against a vacuous comparison: a different seed must NOT match.
add_test(NAME golden_seed43_differs
  COMMAND ${CMAKE_COMMAND}
    -DSIM=$<TARGET_FILE:wormCP2Mu2dAFMplaq>
    -DGOLDEN_DIR=${GOLDEN_DIR}
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/golden_seed43
    -DSEED=43 -DEXPECT_MATCH=OFF
    -P ${CMAKE_CURRENT_SOURCE_DIR}/regression/run_golden.cmake)

set(REGRESSION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/regression)

# A seed is chosen at random when none is given; the recorded run_*.toml must
# replay it exactly.
add_test(NAME replay_recorded_run
  COMMAND ${CMAKE_COMMAND}
    -DSIM=$<TARGET_FILE:wormCP2Mu2dAFMplaq>
    -DCONFIG=${REGRESSION_DIR}/random_seed.toml
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/replay
    -P ${REGRESSION_DIR}/run_provenance.cmake)

# Several runs in one process must match separate processes (no state leaks).
add_test(NAME multi_run_matches_separate_runs
  COMMAND ${CMAKE_COMMAND}
    -DSIM=$<TARGET_FILE:wormCP2Mu2dAFMplaq>
    -DCONFIG=${REGRESSION_DIR}/two_runs.toml
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/multi
    -P ${REGRESSION_DIR}/run_multi.cmake)

# Bad input must be rejected with a clear message (and not crash or run).
function(add_bad_config_test name regex)
  add_test(NAME bad_config_${name}
    COMMAND wormCP2Mu2dAFMplaq ${REGRESSION_DIR}/bad_configs/${name}.toml)
  set_tests_properties(bad_config_${name} PROPERTIES PASS_REGULAR_EXPRESSION "error: .*${regex}")
endfunction()

add_bad_config_test(dt_mismatch "no row for dt=0.1")
add_bad_config_test(missing_probs "cannot open probability file")
add_bad_config_test(short_row "expected 51 values, found 4")
add_bad_config_test(unknown_key "unknown key 'mesurements'")
add_bad_config_test(bad_bc "bc must be")
add_bad_config_test(seed_range "seed must be between")
add_bad_config_test(lt_not_multiple_of_4 "Lt must be a multiple of 4")
add_bad_config_test(inline_dt_mismatch "inline probabilities are for dt=0.1")
add_bad_config_test(inline_missing_key "missing 'p48'")
add_bad_config_test(no_probabilities "no probabilities: set 'probabilities'")
add_bad_config_test(beta_precision "beta may have at most 2 decimal places")
add_bad_config_test(mu_precision "mu3 and mu8 may have at most 3 decimal places")
add_bad_config_test(repeats_zero "repeats must be at least 1")
add_bad_config_test(repeats_same_seed "runs 'r_rep1' and 'r_rep2' have the same parameters and seed 5")
add_bad_config_test(same_seed "runs 'a' and 'b' have the same parameters and seed 5")

add_test(NAME bad_option_jobs
  COMMAND wormCP2Mu2dAFMplaq ${REGRESSION_DIR}/two_runs.toml --jobs 0)
set_tests_properties(bad_option_jobs PROPERTIES PASS_REGULAR_EXPRESSION "error: --jobs must be")

# repeats = N expands into NAME_rep1 .. NAME_repN; --only NAME selects all of them.
add_test(NAME repeats_expand
  COMMAND wormCP2Mu2dAFMplaq ${REGRESSION_DIR}/repeats.toml --dry-run)
set_tests_properties(repeats_expand PROPERTIES
  PASS_REGULAR_EXPRESSION "run 'rep_rep1'.*run 'rep_rep2'.*run 'rep_rep3'.*run 'other'")
add_test(NAME repeats_only_group
  COMMAND wormCP2Mu2dAFMplaq ${REGRESSION_DIR}/repeats.toml --only rep --dry-run)
set_tests_properties(repeats_only_group PROPERTIES
  PASS_REGULAR_EXPRESSION "run 'rep_rep1'.*run 'rep_rep3'" FAIL_REGULAR_EXPRESSION "run 'other'")

# beta and lattice may come from [defaults]; a run can override them.
add_test(NAME defaults_beta_lattice
  COMMAND wormCP2Mu2dAFMplaq ${REGRESSION_DIR}/scan_defaults.toml --dry-run)
set_tests_properties(defaults_beta_lattice PROPERTIES
  PASS_REGULAR_EXPRESSION "run 'from_defaults': bc=closed beta=2.5 lattice=4x4x200.*run 'overridden': bc=closed beta=5 lattice=8x4x400")

# With no 'probabilities' key the nearest probabilities.txt in a parent folder is used.
add_test(NAME master_table_found_in_parent
  COMMAND wormCP2Mu2dAFMplaq ${REGRESSION_DIR}/ancestor/sub/found.toml --dry-run)
set_tests_properties(master_table_found_in_parent PROPERTIES
  PASS_REGULAR_EXPRESSION "probabilities: .*ancestor/probabilities.txt line")

# Physics: tiny lattices (10^6 measurements, a few seconds each) compared with exact results by
# `cp2-analysis check-exact`, so they need uv and the analysis package (uv sync --project analysis).
# Without uv they are listed but disabled. Labelled 'exact': skip them with `ctest -LE exact`, run
# only them with `ctest -L exact`.
set(EXACT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/exact)
find_program(UV_EXECUTABLE uv)
if(NOT UV_EXECUTABLE)
  message(STATUS "uv not found: the exact_* tests are disabled")
endif()
function(add_exact_test name)
  add_test(NAME exact_${name}
    COMMAND ${CMAKE_COMMAND}
      -DSIM=$<TARGET_FILE:wormCP2Mu2dAFMplaq>
      -DUV=${UV_EXECUTABLE}
      -DANALYSIS=${PROJECT_SOURCE_DIR}/analysis
      -DCONFIG=${EXACT_DIR}/exact_runs.toml
      -DCSV=${EXACT_DIR}/exactCP2_results.csv
      -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/exact_${name}
      -DNAME=${name}
      -P ${EXACT_DIR}/run_exact.cmake)
  set_tests_properties(exact_${name} PROPERTIES LABELS exact)
  if(NOT UV_EXECUTABLE)
    set_tests_properties(exact_${name} PROPERTIES DISABLED TRUE)
  endif()
endfunction()

# One per [[run]] in exact/exact_runs.toml.
foreach(name 2x2_closed 2x2_yopen 2x3_yopen 2x3_yopen_mu3_0.03 2x3_yopen_mu3_0.1 2x3_yopen_mu3_0.2)
  add_exact_test(${name})
endforeach()
