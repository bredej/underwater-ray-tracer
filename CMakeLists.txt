cmake_minimum_required(VERSION 3.26)

project(underwater-ray-tracer LANGUAGES CXX)

option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_TESTING "Build tests" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(UnderwaterRayTracer INTERFACE)
add_library(urt::urt ALIAS UnderwaterRayTracer)

target_sources(UnderwaterRayTracer
    INTERFACE
        FILE_SET HEADERS
        BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include
        FILES
            include/underwater-ray-tracer/ray_tracer.hpp
)

get_target_property(HeaderList UnderwaterRayTracer HEADER_SET)
message(STATUS "Header files: ${HeaderList}")

include(GNUInstallDirs)

install(TARGETS UnderwaterRayTracer
    EXPORT UnderwaterRayTracerTargets
    FILE_SET HEADERS
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT UnderwaterRayTracerTargets
  FILE UnderwaterRayTracerConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/UnderwaterRayTracer
)

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if (BUILD_TESTING)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    enable_testing()
    add_subdirectory(test)
endif()
