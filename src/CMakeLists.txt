include_directories(../include)

set(EXEC main bench)

foreach(X IN LISTS EXEC)
    add_executable(${X})
    target_sources(${X} PRIVATE
        ${X}.cpp
        system.cpp
    )

    if(OPT)
        target_sources(${X} PRIVATE odes_opt.cpp)
    else()
        target_sources(${X} PRIVATE odes.cpp)
    endif()

    target_link_libraries(${X} PRIVATE mdspan::mdspan OpenMP::OpenMP_CXX)
    target_compile_features(${X} PRIVATE cxx_std_23)

    target_compile_options(${X} PRIVATE -g -march=native -ftree-vectorize -funroll-loops)
endforeach()