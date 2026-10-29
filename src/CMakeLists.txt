add_executable(wormCP2Mu2dAFMplaq worm/main.cpp worm/runner.cpp worm/lattice.cpp worm/config.cpp worm/probabilities.cpp)
target_link_libraries(wormCP2Mu2dAFMplaq PRIVATE ranlxd tomlplusplus::tomlplusplus)
# Recorded in resolved.toml. Taken when CMake configures, so re-run cmake to refresh it.
set(CP2WORM_GIT_DESCRIBE "unknown")
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --always --dirty
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                  OUTPUT_VARIABLE git_describe RESULT_VARIABLE git_rc
                  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(git_rc EQUAL 0 AND git_describe)
    set(CP2WORM_GIT_DESCRIBE "${git_describe}")
  endif()
endif()
target_compile_definitions(wormCP2Mu2dAFMplaq PRIVATE CP2WORM_GIT_DESCRIBE="${CP2WORM_GIT_DESCRIBE}")

foreach(tool analyseWorm getWormXi getXIutvNFiles)
  add_executable(${tool} tools/${tool}.cpp)
  target_include_directories(${tool} PRIVATE ${PROJECT_SOURCE_DIR}/include)
endforeach()

if(CP2WORM_WARNINGS)
  foreach(tgt wormCP2Mu2dAFMplaq analyseWorm getWormXi getXIutvNFiles)
    target_compile_options(${tgt} PRIVATE -Wall -Wextra)
  endforeach()
endif()
