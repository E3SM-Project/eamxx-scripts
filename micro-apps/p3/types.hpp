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

template <typename Real>
using kokkos_2d_t = Kokkos::View<Real**>;

template <typename Real>
using kokkos_1d_t = Kokkos::View<Real*>;

template <typename Real>
using vector_2d_t = std::vector<std::vector<Real> >;

#endif
