cmake_minimum_required(VERSION 3.25)
project(ODE-tests)
enable_language(CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)
FetchContent_Declare(
    mdspan
    GIT_REPOSITORY https://github.com/kokkos/mdspan
    GIT_TAG stable
)
FetchContent_MakeAvailable(mdspan)

add_subdirectory(src)
target_compile_features(main PRIVATE cxx_std_23)