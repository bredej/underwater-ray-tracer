cmake_minimum_required(VERSION 3.13)

project(urt-test LANGUAGES CXX)


add_executable(urt-test
	urt_test.cpp
)

target_link_libraries(urt-test PRIVATE
    underwater-ray-tracer
    ${CONAN_TARGETS}
)

target_compile_features(urt-test PRIVATE cxx_std_17)

add_test(NAME urt-test COMMAND urt-test)
