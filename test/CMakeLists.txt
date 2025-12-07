cmake_minimum_required(VERSION 3.26)

project(urt-test LANGUAGES CXX)

add_executable(urt-test
    urt_test.cpp
)

target_link_libraries(urt-test PRIVATE
    urt::urt
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(urt-test)
