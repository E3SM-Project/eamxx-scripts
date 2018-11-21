#ifndef MICRO_SED_P3_FUNCTIONS_HPP
#define MICRO_SED_P3_FUNCTIONS_HPP

#include "types.hpp"
#include "scream_pack.hpp"
#include "p3_common.hpp"

namespace p3 {
namespace micro_sed {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for p3.
 * TODO: Need a bit more documentation in this file.
 */

template <typename ScalarT, typename DeviceT>
struct Functions
{

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S> using BigPack = scream::pack::BigPack<S>;
  template <typename S> using SmallPack = scream::pack::SmallPack<S>;
  using IntSmallPack = scream::pack::IntSmallPack;

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  template <typename S>
  using Mask = scream::pack::Mask<BigPack<S>::n>;

  template <typename S>
  using SmallMask = scream::pack::Mask<SmallPack<S>::n>;

  using KT = KokkosTypes<Device>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  using view_1d_table = typename KT::template view_1d_table<Scalar, 150>;
  using view_2d_table = typename KT::template view_2d_table<Scalar, 300, 10>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  using MemberType = typename KT::MemberType;

  struct Table3 {
    IntSmallPack dumii, dumjj;
    SmallPack<Scalar> rdumii, rdumjj, inv_dum3;
  };

  //
  // --------- Functions ---------
  //

  KOKKOS_FUNCTION
  static void lookup(const SmallMask<Scalar>& qr_gt_small, Table3& t,
                     const SmallPack<Scalar>& mu_r, const SmallPack<Scalar>& lamr);

  KOKKOS_FUNCTION
  static Spack apply_table(const SmallMask<Scalar>& qr_gt_small, const view_2d_table& table,
                           const Table3& t);

  // Calculate the step in the region [k_bot, k_top].
  template <Int kdir, int nfield>
  KOKKOS_FUNCTION
  static void calc_first_order_upwind_step(
    const Unmanaged<view_1d<const Spack> >& rho,
    const Unmanaged<view_1d<const Spack> >& inv_rho,
    const Unmanaged<view_1d<const Spack> >& inv_dzq,
    const MemberType& team,
    const Int& nk, const Int& k_bot, const Int& k_top, const Scalar& dt_sub,
    const view_1d_ptr_array<Spack, nfield>& flux,
    const view_1d_ptr_array<Spack, nfield>& V,
    const view_1d_ptr_array<Spack, nfield>& r);

  template <int nfield>
  KOKKOS_FUNCTION
  static void calc_first_order_upwind_step(
    const Unmanaged<view_1d<const Spack> >& rho,
    const Unmanaged<view_1d<const Spack> >& inv_rho,
    const Unmanaged<view_1d<const Spack> >& inv_dzq,
    const MemberType& team,
    const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, const Scalar& dt_sub,
    const view_1d_ptr_array<Spack, nfield>& flux,
    const view_1d_ptr_array<Spack, nfield>& V,
    const view_1d_ptr_array<Spack, nfield>& r);

  // Find the bottom and top of the mixing ratio, e.g., qr. It's worth casing
  // these out in two ways: 1 thread/column vs many, and by kdir.
  KOKKOS_FUNCTION
  static Int find_bottom (
    const MemberType& team,
    const Unmanaged<view_1d<const Scalar> >& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present);

  KOKKOS_FUNCTION
  static Int find_top (
    const MemberType& team,
    const Unmanaged<view_1d<const Scalar> >& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present);
};

} // namespace micro_sed
} // namespace p3

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef KOKKOS_ENABLE_CUDA
# include "p3_functions_table3_impl.hpp"
# include "p3_functions_upwind_impl.hpp"
# include "p3_functions_find_impl.hpp"
#endif

#endif
