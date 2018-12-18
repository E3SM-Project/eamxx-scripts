#ifndef INCLUDE_TYPES
#define INCLUDE_TYPES

#include "micro_kokkos.hpp"

#include <vector>

/*
 * Header contains globally useful types for Scream.
 * The global Int and Real types are defined here along
 * with a type dictionary for accessing commonly-used
 * Kokkos types.
 */

typedef int Int;
#ifdef DOUBLE_PRECISION
typedef double Real;
#else
typedef float Real;
#endif

#if defined KOKKOS_COMPILER_GNU
// See https://github.com/kokkos/kokkos-kernels/issues/129
# define ConstExceptGnu
#else
# define ConstExceptGnu const
#endif

// Turn a View's MemoryTraits (traits::memory_traits) into the equivalent
// unsigned int mask. This is an implementation detail for Unmanaged; see next.
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

// The default device we use.
using DefaultDevice = Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

// Struct for getting useful Kokkos types based on the device
template <typename D=DefaultDevice>
struct KokkosTypes
{
  using Device = D;
  using Layout = Kokkos::LayoutRight;
  using MemSpace = typename Device::memory_space;
  using ExeSpace = typename Device::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExeSpace>;
  using MemberType = typename TeamPolicy::member_type;

  template <typename DataType>
  using view = Kokkos::View<DataType, Layout, Device>;

  template <typename Scalar>
  using view_3d = view<Scalar***>;

  template <typename Scalar>
  using view_2d = view<Scalar**>;

  template <typename Scalar>
  using view_1d = view<Scalar*>;

  template <typename Scalar, int X, int Y>
  using view_2d_table = view<const Scalar[X][Y]>;

  template <typename Scalar, int X>
  using view_1d_table = view<const Scalar[X]>;

  // Our workspace implementation makes this a useful type
  template <typename Scalar, int N>
  using view_1d_ptr_array = Kokkos::Array<Unmanaged<view_1d<Scalar> >*, N>;

  template <typename Scalar, int N>
  using view_1d_ptr_carray = Kokkos::Array<const Unmanaged<view_1d<Scalar> >*, N>;
};

template <typename Scalar>
using vector_2d_t = std::vector<std::vector<Scalar> >;

#endif
