cmake_minimum_required(VERSION 3.10)
project(NumericPlusPlus)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

# Include directories
include_directories(include)

# Core C Library
add_library(core_math src/core/math_utils.c)

# Test Executable (C)
add_executable(test_core tests/test_core.c)
target_link_libraries(test_core core_math)
