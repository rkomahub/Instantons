cmake_minimum_required(VERSION 3.14)
project(MonteCarloInstantons)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)

add_executable(montecarlo ${SOURCES})

target_include_directories(montecarlo PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

find_package(Doxygen)

if(DOXYGEN_FOUND)
  set(DOXYGEN_IN  ${CMAKE_SOURCE_DIR}/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

  add_custom_target(docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
  )
else()
  message(STATUS "Doxygen not found, 'docs' target will not be available")
endif()

enable_testing()

add_executable(test_action
    tests/test_01_action.cpp
    src/core/lattice.cpp
    src/core/metropolis.cpp
    src/core/potential.cpp
)

target_include_directories(test_action PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

add_test(NAME action_tests COMMAND test_action)

add_executable(test_potential_bc
    tests/test_02_potential_and_bc.cpp
    src/core/potential.cpp
)

target_include_directories(test_potential_bc PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

add_test(NAME potential_bc_tests COMMAND test_potential_bc)

add_executable(test_cooling_rule
    tests/test_03_cooling_rule.cpp
    src/core/lattice.cpp
    src/core/metropolis.cpp
    src/core/potential.cpp
)

target_include_directories(test_cooling_rule PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

add_test(NAME cooling_rule_tests COMMAND test_cooling_rule)