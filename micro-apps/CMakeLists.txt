cmake_minimum_required(VERSION 3.3)
cmake_policy(SET CMP0057 NEW)

project(SCREAM_MICRO_APPS CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

function (prc var)
  message ("${var}: ${${var}}")
endfunction ()

include(CTest)

# Scream configuration options
set(SCREAM_DOUBLE_PRECISION TRUE CACHE LOGICAL "Set to double precision (default True)")
set(SCREAM_MIMIC_GPU FALSE CACHE LOGICAL "Mimic GPU to correctness test inter-column parallelism on non-GPU platform")
set(SCREAM_ENABLE_FPE FALSE CACHE LOGICAL "Enable floating point error exception")
set(CHUNKSIZE 0 CACHE INTEGER "Chunk size for impls that want it")
set(SCREAM_PACK_SIZE 16 CACHE INTEGER
  "The number of scalars in a scream::pack::Pack and Mask. Larger packs have good performance on conditional-free loops due to improved caching.")
set(SCREAM_SMALL_PACK_SIZE 16 CACHE INTEGER
  "The number of scalars in a scream::pack::SmallPack and SmallMask. Smaller packs can have better performance in loops with conditionals since more of the packs will have masks with uniform value.")

# Check for valid pack sizes
math(EXPR PACK_MODULO "${SCREAM_PACK_SIZE} % ${SCREAM_SMALL_PACK_SIZE}")
if ( (SCREAM_SMALL_PACK_SIZE GREATER SCREAM_PACK_SIZE) OR (NOT PACK_MODULO EQUAL 0) )
  message(FATAL_ERROR "Invalid small pack size. Needs to be <= pack_size and be a factor of it")
endif()

if (Kokkos_DIR)
  include (${Kokkos_DIR}/kokkos_generated_settings.cmake)
  set (Kokkos_INCLUDE_DIR ${Kokkos_DIR}/include)
  set (Kokkos_LIBRARY_DIR ${Kokkos_DIR}/lib)
  set (Kokkos_LINK_FLAGS ${KOKKOS_LINK_FLAGS})
else ()
  message (FATAL_ERROR "SCREAM_MICRO_APPS requires Kokkos_DIR; this is the base directory of the Kokkos installation.")
endif ()

if ("${KOKKOS_GMAKE_ARCH}" STREQUAL "SKX")
  set(SCREAM_SKX_NO_PACK_DEFAULT 1)
else()
  set(SCREAM_SKX_NO_PACK_DEFAULT ${SCREAM_PACK_SIZE})
endif()

set(SCREAM_SKX_NO_PACK ${SCREAM_SKX_NO_PACK_DEFAULT} CACHE INTEGER
  "Set skx-no-pack to this value. Expect this to usually be 1, but you can set it to something else to restore packs on SKX for testing.")

string (REPLACE ";" " " KOKKOS_CXXFLAGS_STR "${KOKKOS_CXXFLAGS}")
set(KOKKOS_LDFLAGS_STR "")
list(LENGTH KOKKOS_LDFLAGS THE_LEN)
list(LENGTH KOKKOS_CXXFLAGS THE_LEN)
foreach(LDFLAG ${KOKKOS_LDFLAGS})
  if (${LDFLAG} IN_LIST KOKKOS_CXXFLAGS)
    message("Skipping repeat flag ${LDFLAG}")
  else()
    set(KOKKOS_LDFLAGS_STR "${KOKKOS_LDFLAGS_STR} ${LDFLAG}")
  endif()
endforeach()

set (SCREAM_MICRO_APPS_COMPILE_FLAGS ${KOKKOS_CXXFLAGS_STR})
set (SCREAM_MICRO_APPS_LINK_FLAGS ${KOKKOS_LDFLAGS_STR})
set (SCREAM_MICRO_APPS_INCLUDES ${Kokkos_INCLUDE_DIR})
set (SCREAM_MICRO_APPS_LIBRARIES ${KOKKOS_LIBS})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SCREAM_MICRO_APPS_COMPILE_FLAGS} -Wall -fopenmp")

if (${SCREAM_DOUBLE_PRECISION})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DDOUBLE_PRECISION")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOUBLE_PRECISION")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64 ")
  else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 ")
  endif()
endif()

string(FIND "${KOKKOS_GMAKE_DEVICES}" "Cuda" cuda_str_pos)
set(CUDA_BUILD FALSE)
if (${cuda_str_pos} GREATER -1)
  set(CUDA_BUILD TRUE)
endif ()

prc(KOKKOS_GMAKE_DEVICES)
prc(DOUBLE_PRECISION)
prc(CHUNKSIZE)
prc(CUDA_BUILD)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict")
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")
endif()

if ("${KOKKOS_GMAKE_ARCH}" STREQUAL "KNL")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xMIC-AVX512")
elseif ("${KOKKOS_GMAKE_ARCH}" STREQUAL "SKX")
  # Testing shows Fortran does better on SKX with AVX2 than with AVX512.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xCORE-AVX2")
endif()

STRING(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_ci)
if (CMAKE_BUILD_TYPE_ci STREQUAL "debug")
  if (${SCREAM_PACK_SIZE} EQUAL 1)
    set(SCREAM_ENABLE_FPE TRUE)
  endif ()
  set(SCREAM_MIMIC_GPU TRUE)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message("Disabling vectorization in this DEBUG build.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-vec -fp-model strict")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-vec -fp-model strict")
  endif ()
endif ()

if (${SCREAM_ENABLE_FPE})
  message("Enabling FPE.")
  add_definitions(-DFPE)
endif()

if (${SCREAM_MIMIC_GPU})
  if (CUDA_BUILD)
    message("On GPU, so ignore MIMIC_GPU.")
  else ()
    message("Mimic GPU to correctness test inter-column parallelism on non-GPU platform.")
    add_definitions(-DMIMIC_GPU)
  endif ()
endif()

if (CUDA_BUILD)
  set(SCREAM_PACK_SIZE 1)
  set(SCREAM_SMALL_PACK_SIZE 1)
  set(SCREAM_SKX_NO_PACK 1)
endif()

add_definitions(-DCHUNKSIZE=${CHUNKSIZE})
add_definitions(-DSCREAM_PACKN=${SCREAM_PACK_SIZE})
add_definitions(-DSCREAM_SMALL_PACKN=${SCREAM_SMALL_PACK_SIZE})
add_definitions(-DSCREAM_SKX_NO_PACKN=${SCREAM_SKX_NO_PACK})

enable_testing()

add_subdirectory(nano)
add_subdirectory(share)
add_subdirectory(p3)
add_subdirectory(micro-sed)
add_subdirectory(lin-interp)
