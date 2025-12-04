cmake_minimum_required(VERSION 3.10)
project(NumericPlusPlus LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Executable
add_executable(${PROJECT_NAME} main.cpp)
