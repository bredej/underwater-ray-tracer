#!/bin/bash

# To build this project do
# sudo chmod setup.sh +x
# ./setup.sh

# Exit immediately if a command exits with a non-zero status.
set -e

rm -rf build

# Configure CMake with examples and tests enabled
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON -DBUILD_TESTING=ON -B build

# Build the project
cmake --build build

# Run tests
ctest --test-dir build
