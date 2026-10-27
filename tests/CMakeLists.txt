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
