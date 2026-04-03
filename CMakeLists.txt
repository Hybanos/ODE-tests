cmake_minimum_required(VERSION 3.25)
project(ODE-tests)
enable_language(CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(src)
target_compile_features(main PRIVATE cxx_std_23)