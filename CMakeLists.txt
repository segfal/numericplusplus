cmake_minimum_required(VERSION 3.10)
project(NumericPlusPlus LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Headers)

# Executable
add_executable(${PROJECT_NAME} main.cpp)

# Scripts directory
set(SCRIPTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/scripts)

# Custom targets for automation
add_custom_target(clean-script
    COMMAND ${SCRIPTS_DIR}/clean.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running clean script"
)

add_custom_target(build-script
    COMMAND ${SCRIPTS_DIR}/build.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running build script"
)

add_custom_target(run-script
    COMMAND ${SCRIPTS_DIR}/run.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running run script"
)
