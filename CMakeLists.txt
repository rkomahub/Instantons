cmake_minimum_required(VERSION 3.14)
project(MonteCarloInstantons)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Source files
file(GLOB SOURCES CONFIGURE_DEPENDS src/*.cpp)

add_executable(montecarlo ${SOURCES})

# Optional: treat warnings as errors
# target_compile_options(montecarlo PRIVATE -Wall -Wextra -Werror)