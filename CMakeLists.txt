cmake_minimum_required(VERSION 3.19)
project(cp2worm LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(CP2WORM_WARNINGS "Enable extra compiler warnings" OFF)

# Boost is header-only here (format, multi_array). Point CMake at a local
# install with -DBoost_ROOT=... (see CMakeUserPresets.json in the README).
# An unpacked Boost source tree has no BoostConfig.cmake, so keep CMake's FindBoost module
# (deprecated from CMake 3.30, policy CMP0167). Goes away once Boost is no longer needed.
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()
find_package(Boost 1.74 REQUIRED)

add_subdirectory(third_party/ranlxd)
add_subdirectory(third_party/tomlplusplus)
add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
