add_executable(main)

target_sources(main PRIVATE
    main.cpp
    
    system.cpp
    methods.cpp
)

include_directories(../include)
target_link_libraries(main PRIVATE mdspan::mdspan)