cmake_minimum_required(VERSION 3.13)

project(plot LANGUAGES CXX)


add_executable(plot
	plot.cpp
)

target_link_libraries(plot
    underwater-ray-tracer
    ${CONAN_TARGETS}
)

target_compile_features(plot PRIVATE cxx_std_17)

# Suppress compiler warnings from gnuplot-iostream.h
target_compile_definitions(plot PRIVATE $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>)
