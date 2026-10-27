# An unseeded run records its seed; replaying the recorded file must reproduce
# the output byte for byte. Inputs: SIM, CONFIG, WORK_DIR.
file(REMOVE_RECURSE ${WORK_DIR})
foreach(d first replay)
  file(MAKE_DIRECTORY ${WORK_DIR}/${d})
endforeach()

execute_process(COMMAND ${SIM} ${CONFIG} --output-dir ${WORK_DIR}/first
                RESULT_VARIABLE rc OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "first run failed (${rc})")
endif()

file(GLOB record ${WORK_DIR}/first/run_*.toml)
list(LENGTH record n)
if(NOT n EQUAL 1)
  message(FATAL_ERROR "expected exactly one run_*.toml record, found ${n}")
endif()
file(READ ${record} text)
if(NOT text MATCHES "seed = [0-9]+")
  message(FATAL_ERROR "record has no seed:\n${text}")
endif()

execute_process(COMMAND ${SIM} ${record} --output-dir ${WORK_DIR}/replay
                RESULT_VARIABLE rc OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "replay failed (${rc})")
endif()

file(GLOB dats RELATIVE ${WORK_DIR}/first ${WORK_DIR}/first/*.dat)
list(LENGTH dats n)
if(NOT n EQUAL 2)
  message(FATAL_ERROR "expected C_ and M_ files, found ${n}")
endif()
foreach(f ${dats})
  execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${WORK_DIR}/first/${f} ${WORK_DIR}/replay/${f}
                  RESULT_VARIABLE differs)
  if(differs)
    message(FATAL_ERROR "${f} differs between the run and its replay")
  endif()
endforeach()
