# Runs the simulator with a fixed seed and compares its .dat output with the
# golden files. Inputs: SIM, GOLDEN_DIR, WORK_DIR, SEED, EXPECT_MATCH.
#
# To regenerate the golden after an intended change to the physics:
#   wormCP2Mu2dAFMplaq tests/regression/golden/golden.toml --output-dir tests/regression/golden

file(REMOVE_RECURSE ${WORK_DIR})
file(MAKE_DIRECTORY ${WORK_DIR})

execute_process(
  COMMAND ${SIM} ${GOLDEN_DIR}/golden.toml --output-dir ${WORK_DIR} --seed ${SEED}
  RESULT_VARIABLE rc
  OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "simulator exited with ${rc}")
endif()

set(suffix "wormCP2AFMmu2dplaq_BC_closed_nMeas_1000_beta_2.50_XYT_4_4_200_mu3_0.000_mu8_0.000_nSim_1.dat")
set(all_match TRUE)
foreach(kind C M)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files
      ${GOLDEN_DIR}/${kind}_${suffix} ${WORK_DIR}/${kind}_${suffix}
    RESULT_VARIABLE differs
    OUTPUT_QUIET ERROR_QUIET)
  if(differs)
    set(all_match FALSE)
    message(STATUS "${kind}_*.dat differs from golden (seed ${SEED})")
  endif()
endforeach()

if(EXPECT_MATCH AND NOT all_match)
  message(FATAL_ERROR "output does not match golden")
elseif(NOT EXPECT_MATCH AND all_match)
  message(FATAL_ERROR "output matched golden despite a different seed; the comparison is vacuous")
endif()
