#ifndef P3_FINAL_HPP
#define P3_FINAL_HPP

#include "util.hpp"
#include "initial_conditions.hpp"
#include "micro_kokkos.hpp"
#include "p3_common.hpp"
#include "scream_pack.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace p3 {
namespace micro_sed {

using scream::pack::IntPack;
using scream::pack::IntSmallPack;
using scream::pack::smallize;
using scream::pack::scalarize;
using scream::pack::BigPack;
using scream::pack::SmallPack;

template <typename Scalar>
using Mask = scream::pack::Mask<BigPack<Scalar>::n>;

template <typename Scalar>
using SmallMask = scream::pack::Mask<SmallPack<Scalar>::n>;

template <typename Scalar>
struct Table3 {
  IntSmallPack dumii, dumjj;
  SmallPack<Scalar> rdumii, rdumjj, inv_dum3;
};

template <typename Scalar>
KOKKOS_INLINE_FUNCTION
void find_lookupTable_indices_3_kokkos (
  const SmallMask<Scalar>& qr_gt_small, Table3<Scalar>& t, const SmallPack<Scalar>& mu_r, const SmallPack<Scalar>& lamr);

template <typename Scalar, typename D>
struct MicroSedFuncFinalKokkos
{
  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  template <typename S>
  using view_1d = typename KokkosTypes<D>::template view_1d<S>;
  template <typename S>
  using view_2d = typename KokkosTypes<D>::template view_2d<S>;

  using view_1d_table = typename KokkosTypes<D>::template view_1d_table<Scalar, 150>;
  using view_2d_table = typename KokkosTypes<D>::template view_2d_table<Scalar, 300, 10>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KokkosTypes<D>::template view_1d_ptr_carray<S, N>;

  using ExeSpace    = typename KokkosTypes<D>::ExeSpace;
  using MemberType  = typename KokkosTypes<D>::MemberType;
  using TeamPolicy  = typename KokkosTypes<D>::TeamPolicy;

public:
  MicroSedFuncFinalKokkos(int num_horz_, int num_vert_);

  int get_num_vert() const;

  static std::string custom_msg();

  KOKKOS_INLINE_FUNCTION
  static Spack apply_table (
    const SmallMask<Scalar>& qr_gt_small, const view_2d_table& table, const Table3<Scalar>& t);

  // Computes and returns rain size distribution parameters
  KOKKOS_INLINE_FUNCTION
  static void get_rain_dsd2_kokkos (
    const view_1d_table& mu_r_table,
    const SmallMask<Scalar>& qr_gt_small, const Spack& qr, Spack& nr, Spack& mu_r,
    Spack& rdumii, IntSmallPack& dumii, Spack& lamr,
    Spack& cdistr, Spack& logn0r);

  //TODO Unit test.
  // Calculate the step in the region [k_bot, k_top].
  template <Int kdir, int nfield>
  KOKKOS_INLINE_FUNCTION
  static void calc_first_order_upwind_step (
    const Unmanaged<view_1d<const Spack> >& rho,
    const Unmanaged<view_1d<const Spack> >& inv_rho,
    const Unmanaged<view_1d<const Spack> >& inv_dzq,
    const MemberType& team,
    const Int& nk, const Int& k_bot, const Int& k_top, const Scalar& dt_sub,
    const view_1d_ptr_array<Spack, nfield>& flux,
    const view_1d_ptr_array<Spack, nfield>& V,
    const view_1d_ptr_array<Spack, nfield>& r);

  template <int nfield>
  KOKKOS_INLINE_FUNCTION
  static void calc_first_order_upwind_step (
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
  KOKKOS_INLINE_FUNCTION
  static Int find_bottom (
    const MemberType& team,
    const Unmanaged<view_1d<const Scalar> >& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present);

  KOKKOS_INLINE_FUNCTION
  static Int find_top (
    const MemberType& team,
    const Unmanaged<view_1d<const Scalar> >& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present);

  static void micro_sed_func (
    const MicroSedFuncFinalKokkos<Scalar, D>& msfk,
    const Int kts, const Int kte, const int its, const int ite, const Scalar dt,
    const view_2d<Pack>& qr, const view_2d<Pack>& nr,
    const view_2d<Pack>& th, const view_2d<Pack>& dzq, const view_2d<Pack>& pres,
    const view_1d<Scalar>& prt_liq);
};

} // namespace p3
} // namespace micro_sed

#endif
