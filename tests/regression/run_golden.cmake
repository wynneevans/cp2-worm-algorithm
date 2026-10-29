# Runs one [[run]] of golden/golden.toml with a fixed seed and compares its M_ and C_ files
# with golden/M_<NAME>.dat and golden/C_<NAME>.dat, or, where golden/<kind>_<NAME>.sha256 exists
# instead (cases too big to store), with that SHA-256 checksum.
# Inputs: SIM, GOLDEN_DIR, WORK_DIR, NAME, FOLDER (its parameter folder), SEED, EXPECT_MATCH.
#
# To regenerate the golden files after an intended change to the physics:
#   wormCP2Mu2dAFMplaq tests/regression/golden/golden.toml --data-dir /tmp/g --seed 42
#   then copy each /tmp/g/<folder>/M_<time>.dat and C_<time>.dat to tests/regression/golden/M_<name>.dat
#   and C_<name>.dat (the run names and folders are listed in tests/CMakeLists.txt); for the cases with
#   a .sha256 file write `shasum -a 256 <file> | cut -d' ' -f1` into golden/<kind>_<name>.sha256 instead

file(REMOVE_RECURSE ${WORK_DIR})
file(MAKE_DIRECTORY ${WORK_DIR})

execute_process(
  COMMAND ${SIM} ${GOLDEN_DIR}/golden.toml --only ${NAME} --data-dir ${WORK_DIR} --seed ${SEED}
  RESULT_VARIABLE rc
  OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "simulator exited with ${rc}")
endif()

set(folder ${WORK_DIR}/${FOLDER})
set(all_match TRUE)
foreach(kind C M)
  file(GLOB out ${folder}/${kind}_*.dat)
  list(LENGTH out n)
  if(NOT n EQUAL 1)
    message(FATAL_ERROR "expected one ${kind}_*.dat in ${folder}, found ${n}")
  endif()
  set(hash_file ${GOLDEN_DIR}/${kind}_${NAME}.sha256)
  if(EXISTS ${hash_file})
    file(STRINGS ${hash_file} expected LIMIT_COUNT 1)
    file(SHA256 ${out} actual)
    if(actual STREQUAL expected)
      set(differs FALSE)
    else()
      set(differs TRUE)
    endif()
  else()
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E compare_files ${GOLDEN_DIR}/${kind}_${NAME}.dat ${out}
      RESULT_VARIABLE differs
      OUTPUT_QUIET ERROR_QUIET)
  endif()
  if(differs)
    set(all_match FALSE)
    message(STATUS "${kind}_ file differs from golden (${NAME}, seed ${SEED})")
  endif()
endforeach()

if(EXPECT_MATCH AND NOT all_match)
  message(FATAL_ERROR "output does not match golden")
elseif(NOT EXPECT_MATCH AND all_match)
  message(FATAL_ERROR "output matched golden despite a different seed; the comparison is vacuous")
endif()
