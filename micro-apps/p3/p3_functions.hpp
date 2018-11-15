#ifndef MICRO_SED_P3_FUNCTIONS_HPP
#define MICRO_SED_P3_FUNCTIONS_HPP

#include "types.hpp"
#include "scream_pack.hpp"
#include "p3_common.hpp"

namespace p3 {
namespace micro_sed {

template <typename ScalarT, typename DeviceT>
struct Functions {
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

public:  
  struct Table3 {
    IntSmallPack dumii, dumjj;
    SmallPack<Scalar> rdumii, rdumjj, inv_dum3;
  };

  KOKKOS_FUNCTION
  static void lookup(const SmallMask<Scalar>& qr_gt_small, Table3& t,
                     const SmallPack<Scalar>& mu_r, const SmallPack<Scalar>& lamr);

  KOKKOS_FUNCTION
  static Spack apply_table(const SmallMask<Scalar>& qr_gt_small, const view_2d_table& table,
                           const Table3& t);

  //TODO Unit test.
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

  //TODO Unit test.
  // Find the bottom and top of the mixing ratio, e.g., qr. It's worth casing
  // these out in two ways: 1 thread/column vs many, and by kdir.
  KOKKOS_INLINE_FUNCTION
  static Int find_bottom (
    const MemberType& team,
    const Unmanaged<view_1d<const Scalar> >& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present)
  {
    log_present = false;
    Int k_xbot = 0;
    if (team.team_size() == 1) {
      for (Int k = kbot; k != ktop + kdir; k += kdir) {
        if (v(k) < small) continue;
        k_xbot = k;
        log_present = true;
        break;
      }
    } else {
      if (kdir == -1) {
        Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team, kbot - ktop + 1), [&] (Int k_, int& lmax) {
            const Int k = ktop + k_;
            if (v(k) >= small && k > lmax)
              lmax = k;
          }, Kokkos::Max<int>(k_xbot));
        log_present = k_xbot >= ktop;
      } else {
        Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team, ktop - kbot + 1), [&] (Int k_, int& lmin) {
            const Int k = kbot + k_;
            if (v(k) >= small && k < lmin)
              lmin = k;
          }, Kokkos::Min<int>(k_xbot));
        log_present = k_xbot <= ktop;
      }
    }
    return k_xbot;
  }

  //TODO Unit test.
  KOKKOS_INLINE_FUNCTION
  static Int find_top (
    const MemberType& team,
    const Unmanaged<view_1d<const Scalar> >& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present)
  {
    log_present = false;
    Int k_xtop = 0;
    if (team.team_size() == 1) {
      for (Int k = ktop; k != kbot - kdir; k -= kdir) {
        if (v(k) < small) continue;
        k_xtop = k;
        log_present = true;
        break;
      }
    } else {
      if (kdir == -1) {
        Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team, kbot - ktop + 1), [&] (Int k_, int& lmin) {
            const Int k = ktop + k_;
            if (v(k) >= small && k < lmin)
              lmin = k;
          }, Kokkos::Min<int>(k_xtop));
        log_present = k_xtop <= kbot;
      } else {
        Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team, ktop - kbot + 1), [&] (Int k_, int& lmax) {
            const Int k = kbot + k_;
            if (v(k) >= small && k > lmax)
              lmax = k;
          }, Kokkos::Max<int>(k_xtop));
        log_present = k_xtop >= kbot;
      }
    }
    return k_xtop;
  }

  template <typename T>
  void foo(T& t) {}
};

  template <typename T>
  void bar(T& t) {}


} // namespace micro_sed
} // namespace p3

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
//TODO But optional ETI is not yet implemented for these functions.
//#ifdef KOKKOS_ENABLE_CUDA
//# include "p3_functions_table3.hpp"
//# include "p3_functions_upwind.hpp"
//#endif

#endif
