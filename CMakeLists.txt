# This CMake file is more complicated than it needs to be 
# to generate a Visual Studio solution tuned for the developer.
cmake_minimum_required(VERSION 3.13)

project(underwater-ray-tracer LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_TESTING "Build tests" OFF)

if(EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
elseif(EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo_multi.cmake)
    include(${CMAKE_BINARY_DIR}/conanbuildinfo_multi.cmake)
endif()

conan_basic_setup(TARGETS)

set(HEADER_FILES
	"${CMAKE_CURRENT_SOURCE_DIR}/include/underwater-ray-tracer/ray_tracer.hpp"
)

add_library(underwater-ray-tracer INTERFACE)
add_library(underwater-ray-tracer::underwater-ray-tracer ALIAS underwater-ray-tracer)

target_sources(underwater-ray-tracer INTERFACE "$<BUILD_INTERFACE:${HEADER_FILES}>")
target_include_directories(underwater-ray-tracer
	SYSTEM INTERFACE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
		$<INSTALL_INTERFACE:include/>
)

target_compile_features(underwater-ray-tracer INTERFACE cxx_std_17)

if (MSVC)
	# Create dummy target to show header files in IDE
	add_custom_target(underwater-ray-tracer_
		DEPENDS
			${HEADER_FILES}
		SOURCES
			${HEADER_FILES}
	)
endif()


if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if (BUILD_TESTING)
    add_subdirectory(test)
endif()
