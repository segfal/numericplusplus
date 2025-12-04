#!/bin/bash

# Remove build directory (Out-of-source build)
rm -rf build

# Remove CMake generated files (In-source build artifacts)
rm -rf CMakeFiles
rm -f CMakeCache.txt
rm -f cmake_install.cmake
rm -f Makefile

# Remove generated libraries and executables
rm -f libcore_math.a
rm -f test_core

echo "Cleanup complete."
