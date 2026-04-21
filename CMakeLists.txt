cmake_minimum_required(VERSION 3.26)
project(ODE-tests LANGUAGES CXX )
enable_language(C CXX )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(dop853 STATIC IMPORTED)
set_target_properties(dop853 PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/libdop853.a)

include(FetchContent)
FetchContent_Declare(
    mdspan
    GIT_REPOSITORY https://github.com/kokkos/mdspan
    GIT_TAG stable
)
FetchContent_MakeAvailable(mdspan)
# FetchContent_Declare(
#     Kokkos
#     GIT_REPOSITORY https://github.com/kokkos/kokkos
#     GIT_TAG 5.1.0
# )
# FetchContent_MakeAvailable(Kokkos)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(src)

option(OPT "use optimized code")
option(FPOINT_T "float type")
option(A_TOL "absolute tol")
option(R_TOL "relative tol")

if (OPT)
target_compile_definitions(main PRIVATE OPT=${OPT})
endif()

if (FPOINT_T)
target_compile_definitions(main PRIVATE FPOINT_T=${FPOINT_T})
endif()

if (A_TOL)
target_compile_definitions(main PRIVATE A_TOL=${A_TOL})
endif()

if (R_TOL)
target_compile_definitions(main PRIVATE R_TOL=${R_TOL})
endif()
