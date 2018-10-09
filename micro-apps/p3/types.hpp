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

using ExecSpace = Kokkos::DefaultExecutionSpace;

template <typename Real>
using kokkos_2d_t = Kokkos::View<Real**, Layout, MemSpace>;

template <typename Real>
using kokkos_1d_t = Kokkos::View<Real*, Layout, MemSpace>;

// Short name for views
template <typename DataType, typename... Properties>
using ViewType = Kokkos::View<DataType, Layout, Properties...>;

using MemoryManaged   = Kokkos::MemoryTraits<Kokkos::Restrict>;
using MemoryUnmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Restrict>;

// Managed/Unmanaged view
template <typename DataType, typename... Properties>
using ViewManaged = ViewType<DataType, Properties..., MemoryManaged>;
template <typename DataType, typename... Properties>
using ViewUnmanaged = ViewType<DataType, Properties..., MemoryUnmanaged>;

template <typename Real>
using vector_2d_t = std::vector<std::vector<Real> >;

using team_policy = Kokkos::TeamPolicy<>;
using member_type = team_policy::member_type;
//using thread_policy = Kokkos::TeamThreadRange;
//using vector_policy = Kokkos::ThreadVectorRange;

#endif
