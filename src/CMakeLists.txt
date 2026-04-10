include_directories(../include)

add_executable(main)
target_sources(main PRIVATE
    main.cpp
    system.cpp
    odes.cpp
)
target_link_libraries(main PRIVATE mdspan::mdspan)
target_compile_features(main PRIVATE cxx_std_23)

add_executable(bench)
target_sources(bench PRIVATE
    bench.cpp
    system.cpp
    odes.cpp
)
target_link_libraries(bench PRIVATE mdspan::mdspan)
target_compile_features(bench PRIVATE cxx_std_23)