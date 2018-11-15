#ifndef MICRO_SED_P3_FUNCTIONS_UPWIND_HPP
#define MICRO_SED_P3_FUNCTIONS_UPWIND_HPP

#include "p3_functions.hpp"

namespace p3 {
namespace micro_sed {

template <typename S, typename D>
template <Int kdir, int nfield>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_first_order_upwind_step (
  const Unmanaged<view_1d<const Spack> >& rho,
  const Unmanaged<view_1d<const Spack> >& inv_rho,
  const Unmanaged<view_1d<const Spack> >& inv_dzq,
  const MemberType& team,
  const Int& nk, const Int& k_bot, const Int& k_top, const Scalar& dt_sub,
  const view_1d_ptr_array<Spack, nfield>& flux,
  const view_1d_ptr_array<Spack, nfield>& V,
  const view_1d_ptr_array<Spack, nfield>& r)
{
  Int
    kmin = ( kdir == 1 ? k_bot : k_top)             / Spack::n,
    // Add 1 to make [kmin, kmax). But then the extra term (Spack::n -
    // 1) to determine pack index cancels the +1.
    kmax = ((kdir == 1 ? k_top : k_bot) + Spack::n) / Spack::n;
  const Int k_top_pack = k_top / Spack::n;

  // calculate fluxes
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, kmax - kmin), [&] (Int k_) {
      const Int k = kmin + k_;
      for (int f = 0; f < nfield; ++f)
        (*flux[f])(k) = (*V[f])(k) * (*r[f])(k) * rho(k);
    });
  team.team_barrier();

  Kokkos::single(
    Kokkos::PerTeam(team), [&] () {
      const Int k = k_top_pack;
      if (nk % Spack::n != 0) {
        const auto mask =
          scream::pack::range<IntSmallPack>(k_top_pack*Spack::n) >= nk;
        for (int f = 0; f < nfield; ++f)
          (*flux[f])(k_top_pack).set(mask, 0);
      }
      for (int f = 0; f < nfield; ++f) {
        // compute flux divergence
        const auto flux_pkdir = (kdir == -1) ?
          shift_right(0, (*flux[f])(k)) :
          shift_left (0, (*flux[f])(k));
        const auto fluxdiv = (flux_pkdir - (*flux[f])(k)) * inv_dzq(k);
        // update prognostic variables
        (*r[f])(k) += fluxdiv * dt_sub * inv_rho(k);
      }
    });

  if (kdir == 1)
    --kmax;
  else
    ++kmin;

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, kmax - kmin), [&] (Int k_) {
      const Int k = kmin + k_;
      for (int f = 0; f < nfield; ++f) {
        // compute flux divergence
        const auto flux_pkdir = (kdir == -1) ?
          shift_right((*flux[f])(k+kdir), (*flux[f])(k)) :
          shift_left ((*flux[f])(k+kdir), (*flux[f])(k));
        const auto fluxdiv = (flux_pkdir - (*flux[f])(k)) * inv_dzq(k);
        // update prognostic variables
        (*r[f])(k) += fluxdiv * dt_sub * inv_rho(k);
      }
    });
}

template <typename S, typename D>
template <int nfield>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_first_order_upwind_step (
  const Unmanaged<view_1d<const Spack> >& rho,
  const Unmanaged<view_1d<const Spack> >& inv_rho,
  const Unmanaged<view_1d<const Spack> >& inv_dzq,
  const MemberType& team,
  const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, const Scalar& dt_sub,
  const view_1d_ptr_array<Spack, nfield>& flux,
  const view_1d_ptr_array<Spack, nfield>& V,
  const view_1d_ptr_array<Spack, nfield>& r)
{
  if (kdir == 1)
    calc_first_order_upwind_step< 1, nfield>(
      rho, inv_rho, inv_dzq, team, nk, k_bot, k_top, dt_sub, flux, V, r);
  else
    calc_first_order_upwind_step<-1, nfield>(
      rho, inv_rho, inv_dzq, team, nk, k_bot, k_top, dt_sub, flux, V, r);
}

} // namespace micro_sed
} // namespace p3

#endif
