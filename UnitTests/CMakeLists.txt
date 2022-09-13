cmake_minimum_required(VERSION 3.23)
project(CUniMolSys)
set(CMAKE_CXX_STANDARD 17)

file(GLOB all)

enable_testing()
include(GoogleTest)
add_executable(TestingTest TestingTest.cpp)
target_link_libraries(TestingTest PRIVATE gtest gtest_main)
gtest_discover_tests(TestingTest)