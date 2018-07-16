cmake_minimum_required(VERSION 2.8.11)

project(SCREAM_MICRO_APPS CXX Fortran)

include(CTest)

set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Select build type.")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")

enable_testing()

add_subdirectory(p3)
