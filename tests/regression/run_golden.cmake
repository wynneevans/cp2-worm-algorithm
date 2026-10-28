# Runs the simulator with a fixed seed and compares its M_ and C_ files with the
# golden files. Inputs: SIM, GOLDEN_DIR, WORK_DIR, SEED, EXPECT_MATCH.
#
# To regenerate the golden after an intended change to the physics:
#   wormCP2Mu2dAFMplaq tests/regression/golden/golden.toml --data-dir /tmp/g --seed 42
#   then copy /tmp/g/<folder>/M_<time>.dat and C_<time>.dat to tests/regression/golden/M.dat and C.dat

file(REMOVE_RECURSE ${WORK_DIR})
file(MAKE_DIRECTORY ${WORK_DIR})

execute_process(
  COMMAND ${SIM} ${GOLDEN_DIR}/golden.toml --data-dir ${WORK_DIR} --seed ${SEED}
  RESULT_VARIABLE rc
  OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "simulator exited with ${rc}")
endif()

set(folder ${WORK_DIR}/closed_1000_2.50_4_4_200_mu3_0.000_mu8_0.000)
set(all_match TRUE)
foreach(kind C M)
  file(GLOB out ${folder}/${kind}_*.dat)
  list(LENGTH out n)
  if(NOT n EQUAL 1)
    message(FATAL_ERROR "expected one ${kind}_*.dat in ${folder}, found ${n}")
  endif()
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files ${GOLDEN_DIR}/${kind}.dat ${out}
    RESULT_VARIABLE differs
    OUTPUT_QUIET ERROR_QUIET)
  if(differs)
    set(all_match FALSE)
    message(STATUS "${kind}_ file differs from golden (seed ${SEED})")
  endif()
endforeach()

if(EXPECT_MATCH AND NOT all_match)
  message(FATAL_ERROR "output does not match golden")
elseif(NOT EXPECT_MATCH AND all_match)
  message(FATAL_ERROR "output matched golden despite a different seed; the comparison is vacuous")
endif()
