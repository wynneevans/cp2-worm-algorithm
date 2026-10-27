# Martin Luescher's RANLUX generator (GPL), kept verbatim.
add_library(ranlxd STATIC ranlxd.c)
target_include_directories(ranlxd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
