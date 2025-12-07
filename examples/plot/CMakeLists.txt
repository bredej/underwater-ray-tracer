cmake_minimum_required(VERSION 3.26)

project(plot LANGUAGES CXX)

add_executable(plot
    plot.cpp
)

target_link_libraries(plot PRIVATE
    urt::urt
)
