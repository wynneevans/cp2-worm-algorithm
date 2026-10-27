add_executable(wormCP2Mu2dAFMplaq worm/main.cpp)
target_link_libraries(wormCP2Mu2dAFMplaq PRIVATE Boost::headers ranlxd)
# Boost 1.74's multi_array needs std::unary_function, which libc++ drops in
# C++17. Not needed with Boost >= 1.81; remove once Boost is upgraded.
target_compile_definitions(wormCP2Mu2dAFMplaq PRIVATE
  _LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION)

foreach(tool analyseWorm getWormXi getXIutvNFiles)
  add_executable(${tool} tools/${tool}.cpp)
  target_include_directories(${tool} PRIVATE ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(${tool} PRIVATE Boost::headers)
endforeach()

if(CP2WORM_WARNINGS)
  foreach(tgt wormCP2Mu2dAFMplaq analyseWorm getWormXi getXIutvNFiles)
    target_compile_options(${tgt} PRIVATE -Wall -Wextra)
  endforeach()
endif()
