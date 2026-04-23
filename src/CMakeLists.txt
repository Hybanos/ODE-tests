include_directories(../include)

set(EXEC main bench)

foreach(X IN LISTS EXEC)
    add_executable(${X})
    target_sources(${X} PRIVATE
        ${X}.cpp
        system.cpp
        dop853_ref.cpp
    )

    if(OPT)
        target_sources(${X} PRIVATE odes_opt.cpp)
    else()
        target_sources(${X} PRIVATE odes.cpp)
    endif()

    target_link_libraries(${X} PRIVATE mdspan::mdspan OpenMP::OpenMP_CXX)
    target_link_libraries(${X} PRIVATE dop853)
    target_compile_features(${X} PRIVATE cxx_std_23)

    target_compile_options(${X} PRIVATE -g -march=native -ftree-vectorize -fno-omit-frame-pointer -grecord-gcc-switches)
    # if(OPT)
    #     target_compile_options(${X} PRIVATE -ffast-math)
    # endif()
endforeach()