cmake_minimum_required(VERSION 3.10)
project(CUniMolSys)
set(CMAKE_CXX_STANDARD 20)

if(${MSVC})
    add_compile_options("/MT")
    add_link_options("/MACHINE:X64")
    add_link_options("/LARGEADDRESSAWARE") # without this, MVSC will fail on >2GB memory usage.
    IF(CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options("/MTd")
    else()
        add_compile_options("/MT")
    ENDIF()
endif()

add_compile_options("-Wall")

# Enable GoogleTests
add_subdirectory(googletest)

# Configure the DataFilesPath header "resources.h" from "resources.in"
set(_DATAFILESPATH_ \"${CMAKE_SOURCE_DIR}/Data\")
message("${_DATAFILESPATH_}")
configure_file(resources.in ${CMAKE_SOURCE_DIR}/resources.h)

# Automatically find all .cpp and .h files in current folder
file(GLOB SRCS ${CMAKE_CURRENT_LIST_DIR}/*.cpp)
file(GLOB HEADERS ${CMAKE_CURRENT_LIST_DIR}/*.h)

# exclude main.cpp from LIBSOURCE
set(LIBSOURCE "")
foreach(source IN LISTS SRCS)
    if(NOT source MATCHES "main.cpp")
        set(LIBSOURCE "${LIBSOURCE};${source}")
    endif()
endforeach()
#message(${LIBSOURCE})
#message(${HEADERS})

# Copy *.h to ./include
foreach(header IN LISTS HEADERS)
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${header} ${CMAKE_CURRENT_LIST_DIR}/include/)
endforeach()

set(LIBNAME cunimolsys)
add_library(${LIBNAME} STATIC "${LIBSOURCE}" "${HEADERS}")
set_target_properties(${LIBNAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)

add_executable(main main.cpp)
target_link_libraries(main ${LIBNAME})


# Add All Tests
add_subdirectory(UnitTests)

# Add applications
add_subdirectory(TrajAnalysis)
