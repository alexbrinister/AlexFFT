cmake_minimum_required(VERSION 3.10)

project(dft_program)

# Set the C standard to C17
set(CMAKE_C_STANDARD 17)

# C flags
set(CFLAGS CFLAGS -Wall -Werror -pedantic -Wextra -fPIE -fstack-protector-all)

# Add the executable
add_executable(dft_program src/fft.c)

# Link the math library
target_link_libraries(dft_program m)
