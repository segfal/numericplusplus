cmake_minimum_required(VERSION 3.10)
project(C_Numerical_Methods LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add the new C project subdirectory
add_subdirectory(damped_oscillator)

set(SCRIPTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/scripts)

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

# Update the run-script to target the new executable
add_custom_target(run-script
    COMMAND ${SCRIPTS_DIR}/run.sh damped_oscillator
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running damped_oscillator"
)
