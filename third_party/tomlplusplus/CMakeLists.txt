# toml++ single header, vendored unmodified (see VERSION.txt).
add_library(tomlplusplus INTERFACE)
add_library(tomlplusplus::tomlplusplus ALIAS tomlplusplus)
target_include_directories(tomlplusplus SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
