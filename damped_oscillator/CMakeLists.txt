# damped_oscillator/CMakeLists.txt

add_executable(damped_oscillator main.c oscillator.c)

# Link the standard math library (libm)
target_link_libraries(damped_oscillator PRIVATE m)
