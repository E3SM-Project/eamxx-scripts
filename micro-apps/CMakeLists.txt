cmake_minimum_required(VERSION 2.8.11)

project(SCREAM_MICRO_APPS CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

function (prc var)
  message ("${var}: ${${var}}")
endfunction ()

include(CTest)

set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Select build type.")
set(DOUBLE_PRECISION TRUE CACHE LOGICAL "Set to double precision (default True)")
set(ENABLE_TRACING FALSE CACHE LOGICAL "Enable detailed tracing")
set(ENABLE_FPE FALSE CACHE LOGICAL "Enable floating point error exception")

if (Kokkos_DIR)
  include (${Kokkos_DIR}/kokkos_generated_settings.cmake)
  set (Kokkos_INCLUDE ${Kokkos_DIR}/include)
  set (Kokkos_LIBRARY_DIR ${Kokkos_DIR}/lib)
else ()
  message (FATAL_ERROR "SCREAM_MICRO_APPS requires Kokkos_DIR; this is the base directory of the Kokkos installation.")
endif ()

string (REPLACE ";" " " KOKKOS_CXXFLAGS_STR "${KOKKOS_CXXFLAGS}")
string (REPLACE ";" " " KOKKOS_LDFLAGS_STR "${KOKKOS_LDFLAGS}")

set (SCREAM_MICRO_APPS_COMPILE_FLAGS ${KOKKOS_CXXFLAGS_STR})
set (SCREAM_MICRO_APPS_LINK_FLAGS ${KOKKOS_LDFLAGS_STR})
set (SCREAM_MICRO_APPS_INCLUDES ${Kokkos_INCLUDE})
set (SCREAM_MICRO_APPS_LIBRARIES ${KOKKOS_LIBS})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -cpp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SCREAM_MICRO_APPS_COMPILE_FLAGS} -Wall")

if (${DOUBLE_PRECISION})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -D DOUBLE_PRECISION")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D DOUBLE_PRECISION")
endif()

if (${ENABLE_TRACING})
  add_definitions(-DTRACE)
endif()

if (${ENABLE_FPE})
  add_definitions(-DFPE)
endif()

enable_testing()

add_subdirectory(p3)
