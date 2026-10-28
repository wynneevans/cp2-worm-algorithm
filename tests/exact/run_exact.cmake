# Runs one [[run]] of exact_runs.toml and checks it against the exact results with
# `cp2-analysis check-exact`, so the same estimators as the analysis are tested.
# Inputs: SIM, UV, ANALYSIS, CONFIG, CSV, WORK_DIR, NAME.
file(REMOVE_RECURSE ${WORK_DIR})
execute_process(COMMAND ${SIM} ${CONFIG} --only ${NAME} --data-dir ${WORK_DIR}
                RESULT_VARIABLE rc OUTPUT_QUIET)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "simulator exited with ${rc}")
endif()
execute_process(COMMAND ${UV} run --quiet --project ${ANALYSIS} cp2-analysis check-exact ${WORK_DIR}
                        --typed --exact ${CSV}
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "Monte Carlo disagrees with the exact results (or cp2-analysis failed: ${rc})")
endif()
