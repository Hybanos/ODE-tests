cmake_minimum_required(VERSION 3.26)
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

option(FPOINT_T "float type")
option(A_TOL "absolute tol")
option(R_TOL "relative tol")

if (FPOINT_T)
target_compile_definitions(main PRIVATE FPOINT_T=${FPOINT_T}
)
endif()

if (A_TOL)
target_compile_definitions(main PRIVATE A_TOL=${A_TOL}
)
endif()

if (R_TOL)
target_compile_definitions(main PRIVATE R_TOL=${R_TOL}
)
endif()
