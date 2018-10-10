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

// Turn a View's MemoryTraits (traits::memory_traits) into the equivalent
// unsigned int mask.
template <typename View>
struct MemoryTraitsMask {
  enum : unsigned int {
    value = ((View::traits::memory_traits::RandomAccess ? Kokkos::RandomAccess : 0) |
             (View::traits::memory_traits::Atomic ? Kokkos::Atomic : 0) |
             (View::traits::memory_traits::Restrict ? Kokkos::Restrict : 0) |
             (View::traits::memory_traits::Aligned ? Kokkos::Aligned : 0) |
             (View::traits::memory_traits::Unmanaged ? Kokkos::Unmanaged : 0))
      };
};

// Make the input View Unmanaged, whether or not it already is. One might
// imagine that View::unmanaged_type would provide this.
//   Use: Unmanged<ViewType>
template <typename View>
using Unmanaged =
  // Provide a full View type specification, augmented with Unmanaged.
  Kokkos::View<typename View::traits::scalar_array_type,
               typename View::traits::array_layout,
               typename View::traits::device_type,
               Kokkos::MemoryTraits<
                 // All the current values...
                 MemoryTraitsMask<View>::value |
                 // ... |ed with the one we want, whether or not it's
                 // already there.
                 Kokkos::Unmanaged> >;

template <typename Real>
using vector_2d_t = std::vector<std::vector<Real> >;

using team_policy = Kokkos::TeamPolicy<>;
using member_type = team_policy::member_type;
//using thread_policy = Kokkos::TeamThreadRange;
//using vector_policy = Kokkos::ThreadVectorRange;

#endif
