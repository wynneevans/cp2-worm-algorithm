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
