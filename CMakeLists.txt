cmake_minimum_required(VERSION 3.19)
project(cp2worm LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(CP2WORM_WARNINGS "Enable extra compiler warnings" OFF)

# Boost is header-only here (format, multi_array). Point CMake at a local
# install with -DBoost_ROOT=... (see CMakeUserPresets.json in the README).
find_package(Boost 1.74 REQUIRED)

add_subdirectory(third_party/ranlxd)
add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
