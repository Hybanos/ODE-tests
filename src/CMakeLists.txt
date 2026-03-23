add_executable(main)

target_sources(main PRIVATE
    main.cpp
    
    system.cpp
    methods/euler.cpp
)

include_directories(../include)