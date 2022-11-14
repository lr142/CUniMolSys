project(TrajAnalysis)
set(CMAKE_CXX_STANDARD 20)


function(my_add_exec name)
    file(GLOB sources *.cpp)
    file(GLOB headers *.h)
    add_executable(${name} "${sources}" "${headers}")
    target_link_libraries(${name} PRIVATE cunimolsys)
    target_link_directories(${name} PRIVATE ../lib)
    target_include_directories(${name} PRIVATE ../include)
endfunction()

if(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    add_definitions(-DDEBUGGING)
endif()

my_add_exec(trajana)