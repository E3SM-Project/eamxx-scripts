#ifndef INCLUDE_TYPES
#define INCLUDE_TYPES

#include "micro_kokkos.hpp"

#include <vector>

typedef int Int;
#ifdef DOUBLE_PRECISION
typedef double Real;
#else
typedef float Real;
#endif

using Layout = Kokkos::LayoutRight;

#ifdef KOKKOS_ENABLE_CUDA
using MemSpace = Kokkos::CudaSpace;
#else
using MemSpace = Kokkos::HostSpace;
#endif

using ExecSpace = Kokkos::Cuda;

template <typename Real>
using kokkos_2d_t = Kokkos::View<Real**, Layout, MemSpace>;

template <typename Real>
using kokkos_1d_t = Kokkos::View<Real*, Layout, MemSpace>;

template <typename Real>
using vector_2d_t = std::vector<std::vector<Real> >;

#endif
