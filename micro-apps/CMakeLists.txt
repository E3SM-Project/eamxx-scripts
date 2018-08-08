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
  include (${Kokkos_DIR}/kokkos.cmake)
  set (Kokkos_INCLUDE ${Kokkos_DIR}/include)
  set (Kokkos_LIBRARY_DIR ${Kokkos_DIR}/lib)
else ()
  message (FATAL_ERROR "COMPOSE requires Kokkos_DIR")
endif ()
set (SCREAM_MICRO_APPS_COMPILE_FLAGS "${KOKKOS_CXXFLAGS}")
set (SCREAM_MICRO_APPS_LINK_FLAGS "${KOKKOS_LDFLAGS}")
set (SCREAM_MICRO_APPS_INCLUDES "${Kokkos_INCLUDE}")
set (SCREAM_MICRO_APPS_LIBRARIES "${KOKKOS_LIBS}")

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
