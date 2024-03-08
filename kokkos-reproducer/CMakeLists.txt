# This nano-project is for creating small reproducers for issues involving Kokkos
# Right now, the reproducer is set to grab the e3sm-kokkos-4.2.00 branch from
# e3sm-project/kokkos, and build a small exec that does an MPI bcast using
# device pointers.

cmake_minimum_required(VERSION 3.14)
include (FetchContent)

project (KokkosIssueReproducer CXX C)

FetchContent_Declare(
  Kokkos
  GIT_REPOSITORY https://github.com/e3sm-project/kokkos.git
  GIT_TAG        e3sm-kokkos-4.2.00
)

# Uncomment the kokkos device(s) you need/want
# set (Kokkos_ENABLE_HIP ON)
# set (Kokkos_ENABLE_CUDA ON)
# set (Kokkos_ENABLE_OPENMP ON)
set (Kokkos_ENABLE_SERIAL ON)

# Uncomment the kokkos arch(s) you need/want
# set (Kokkos_ARCH_VOLTA70 ON)
# set (Kokkos_ARCH_AMPERE80 ON)
# set (Kokkos_ARCH_ZEN2 ON)
# set (Kokkos_ARCH_ZEN3 ON)
# set (Kokkos_ARCH_SKX ON)

# MISC Kokkos settings, uncomment what you need
set (Kokkos_ENABLE_DEPRECATED_CODE OFF)
set (Kokkos_ENABLE_CUDA_LAMBDA ON)


FetchContent_MakeAvailable(Kokkos)

find_package(MPI REQUIRED COMPONENTS C)

add_executable (main main.cpp)
target_link_libraries (main PUBLIC Kokkos::kokkos MPI::MPI_C)
