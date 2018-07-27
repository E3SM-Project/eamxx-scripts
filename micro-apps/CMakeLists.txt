cmake_minimum_required(VERSION 2.8.11)

project(SCREAM_MICRO_APPS CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

include(CTest)

set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Select build type.")
set(DOUBLE_PRECISION TRUE CACHE LOGICAL "Set to double precision (default True)")
set(ENABLE_TRACING FALSE CACHE LOGICAL "Enable detailed tracing")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -cpp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

if (${DOUBLE_PRECISION})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -D DOUBLE_PRECISION")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D DOUBLE_PRECISION")
endif()

if (${ENABLE_TRACING})
  add_definitions(-DTRACE)
endif()

enable_testing()

add_subdirectory(p3)
