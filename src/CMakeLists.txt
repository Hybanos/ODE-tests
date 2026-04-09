add_executable(main2)

target_sources(main2 PRIVATE
    main2.cpp
    
    system2.cpp
    odes.cpp
)

include_directories(../include)
target_link_libraries(main2 PRIVATE mdspan::mdspan)