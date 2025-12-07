@echo off

@rem To build this project do
@rem ./setup.bat

rmdir /S /Q build

@rem Configure CMake with examples and tests enabled
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON -DBUILD_TESTING=ON -B build

@rem Build the project
cmake --build build
