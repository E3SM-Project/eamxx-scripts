#ifndef MICRO_SED_PACKNOIWS_KOKKOS_HPP
#define MICRO_SED_PACKNOIWS_KOKKOS_HPP

#include "util.hpp"
#include "initial_conditions.hpp"
#include "micro_kokkos.hpp"
#include "micro_sed_vanilla.hpp"
#include "micro_sed_vanilla_kokkos.hpp"
#include "micro_sed_packnoi_kokkos.hpp"
#include "scream_pack.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace p3 {
namespace micro_sed {

using micro_sed::Globals;

template <typename Real>
struct MicroSedFuncPackNoiWsKokkos {

  static constexpr const char* NAME = "kokkos_packnoiws";

  int num_horz, num_vert, num_pack;

  kokkos_2d_table_t<Real> vn_table, vm_table;
  kokkos_1d_table_t<Real> mu_r_table;

  team_policy policy;
  util::WorkSpace workspace;

  using pack_t = RealPack;

public:
  MicroSedFuncPackNoiWsKokkos(int num_horz_, int num_vert_) :
    num_horz(num_horz_), num_vert(num_vert_),
    num_pack(scream::pack::npack<RealPack>(num_vert_)),
    vn_table("VN_TABLE"), vm_table("VM_TABLE"),
    mu_r_table("MU_R_TABLE"),
    policy(util::ExeSpaceUtils<>::get_default_team_policy(num_horz, num_pack)),
    workspace(num_vert_ * sizeof(Real), num_horz_, policy)
  {
    // initialize on host

    auto mirror_vn_table = Kokkos::create_mirror_view(vn_table);
    auto mirror_vm_table = Kokkos::create_mirror_view(vm_table);
    auto mirror_mu_table = Kokkos::create_mirror_view(mu_r_table);

    for (int i = 0; i < 300; ++i) {
      for (int k = 0; k < 10; ++k) {
        mirror_vn_table(i, k) = Globals<Real>::VN_TABLE[i][k];
        mirror_vm_table(i, k) = Globals<Real>::VM_TABLE[i][k];
      }
    }

    for (int i = 0; i < 150; ++i) {
      mirror_mu_table(i) = Globals<Real>::MU_R_TABLE[i];
    }

    // deep copy to device
    Kokkos::deep_copy(vn_table, mirror_vn_table);
    Kokkos::deep_copy(vm_table, mirror_vm_table);
    Kokkos::deep_copy(mu_r_table, mirror_mu_table);
  }

  int get_num_vert() const { return num_vert; }

  static std::string custom_msg()
  {
    std::ostringstream out;
    out << " packn=" << SCREAM_PACKN << " small_pack_factor=" << SCREAM_SMALL_PACK_FACTOR;
    return out.str();
  }
};

#if 0
//TODO Unit test.
// Calculate the step in the region [k_bot, k_top].
template <Int kdir, int nfield>
KOKKOS_INLINE_FUNCTION
void calc_first_order_upwind_step (
  const Unmanaged<kokkos_1d_t<const RealSmallPack> >& rho,
  const Unmanaged<kokkos_1d_t<const RealSmallPack> >& inv_rho,
  const Unmanaged<kokkos_1d_t<const RealSmallPack> >& inv_dzq,
  const member_type& team,
  const Int& nk, const Int& k_bot, const Int& k_top, const Real& dt_sub,
  const Kokkos::Array<const Unmanaged<kokkos_1d_t<RealSmallPack> >*,nfield>& flux,
  const Kokkos::Array<const Unmanaged<kokkos_1d_t<RealSmallPack> >*,nfield>& V,
  const Kokkos::Array<const Unmanaged<kokkos_1d_t<RealSmallPack> >*,nfield>& r)
{
  Int
    kmin = ( kdir == 1 ? k_bot : k_top)                     / RealSmallPack::n,
    // Add 1 to make [kmin, kmax). But then the extra term (RealSmallPack::n -
    // 1) to determine pack index cancels the +1.
    kmax = ((kdir == 1 ? k_top : k_bot) + RealSmallPack::n) / RealSmallPack::n;
  const Int k_top_pack = k_top / RealSmallPack::n;

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
      if (nk % RealSmallPack::n != 0) {
        const auto mask =
          scream::pack::range<IntSmallPack>(k_top_pack*RealSmallPack::n) >= nk;
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

template <int nfield> KOKKOS_INLINE_FUNCTION
void calc_first_order_upwind_step (
  const Unmanaged<kokkos_1d_t<const RealSmallPack> >& rho,
  const Unmanaged<kokkos_1d_t<const RealSmallPack> >& inv_rho,
  const Unmanaged<kokkos_1d_t<const RealSmallPack> >& inv_dzq,
  const member_type& team,
  const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, const Real& dt_sub,
  const Kokkos::Array<const Unmanaged<kokkos_1d_t<RealSmallPack> >*,nfield>& flux,
  const Kokkos::Array<const Unmanaged<kokkos_1d_t<RealSmallPack> >*,nfield>& V,
  const Kokkos::Array<const Unmanaged<kokkos_1d_t<RealSmallPack> >*,nfield>& r)
{
  if (kdir == 1)
    calc_first_order_upwind_step< 1, nfield>(
      rho, inv_rho, inv_dzq, team, nk, k_bot, k_top, dt_sub, flux, V, r);
  else
    calc_first_order_upwind_step<-1, nfield>(
      rho, inv_rho, inv_dzq, team, nk, k_bot, k_top, dt_sub, flux, V, r);
}
#endif

void micro_sed_func (
  MicroSedFuncPackNoiWsKokkos<Real>& m,
  const Int kts, const Int kte, const int its, const int ite, const Real dt,
  const kokkos_2d_t<RealPack>& qr, const kokkos_2d_t<RealPack>& nr,
  const kokkos_2d_t<RealPack>& th, const kokkos_2d_t<RealPack>& dzq, const kokkos_2d_t<RealPack>& pres,
  const kokkos_1d_t<Real>& prt_liq)
{
  const kokkos_2d_t<RealSmallPack>
    linv_dzq, //= smallize(m.inv_dzq),
    lrho, //= smallize(m.rho),
    linv_rho, //= smallize(m.inv_rho),
    lrhofacr, //= smallize(m.rhofacr),
    lV_qr, //= smallize(m.V_qr),
    lV_nr, //= smallize(m.V_nr),
    lflux_qx, //= smallize(m.flux_qx),
    lflux_nx, //= smallize(m.flux_nx),
    lqr = smallize(qr),
    lnr = smallize(nr),
    lmu_r, //= smallize(m.mu_r),
    llamr; // = smallize(m.lamr);
  const kokkos_2d_t<Real>
    sqr = scalarize(qr),
    sflux_qx; // = scalarize(m.flux_qx);

  // constants
  const Real odt = 1.0 / dt;
  constexpr auto nsmall = Globals<Real>::NSMALL;
  constexpr auto rd = Globals<Real>::RD;
  constexpr auto rd_inv_cp = Globals<Real>::RD * Globals<Real>::INV_CP;
  constexpr auto rhosur = Globals<Real>::RHOSUR;
  constexpr auto qsmall = Globals<Real>::QSMALL;

  // direction of vertical leveling
  const Int kbot = (kts < kte) ? 0 : m.num_vert-1;
  const Int ktop = (kts < kte) ? m.num_vert-1 : 0;
  const Int kdir = (kts < kte) ? 1 : -1;

  // Rain sedimentation:  (adaptive substepping)
  Kokkos::parallel_for(
    "main rain sed loop",
    m.policy,
    KOKKOS_LAMBDA(const member_type& team) {
      const int i = team.league_rank();

      const Unmanaged<kokkos_1d_t<RealSmallPack> >
        olinv_dzq = util::subview(linv_dzq, i),
        olrho = util::subview(lrho, i),
        olinv_rho = util::subview(linv_rho, i),
        olrhofacr = util::subview(lrhofacr, i),
        olV_qr = util::subview(lV_qr, i),
        olV_nr = util::subview(lV_nr, i),
        olflux_qx = util::subview(lflux_qx, i),
        olflux_nx = util::subview(lflux_nx, i),
        olqr = util::subview(lqr, i),
        olnr = util::subview(lnr, i),
        olmu_r = util::subview(lmu_r, i),
        ollamr = util::subview(llamr, i);

      const Unmanaged<kokkos_1d_t<Real> >
        osqr = util::subview(sqr, i),
        osflux_qx = util::subview(sflux_qx, i);

      const Unmanaged<kokkos_1d_t<RealPack> >
        oinv_dzq,// = util::subview(m.inv_dzq, i),
        ot,     //  = util::subview(m.t, i),
        orho,   //  = util::subview(m.rho, i),
        oinv_rho, //= util::subview(m.inv_rho, i),
        orhofacr, //= util::subview(m.rhofacr, i),
        oV_qr, //   = util::subview(m.V_qr, i),
        oV_nr, //   = util::subview(m.V_nr, i),
        odzq     = util::subview(dzq, i),
        oth      = util::subview(th, i),
        opres    = util::subview(pres, i);

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, m.num_pack), [&] (Int k) {
          oinv_dzq(k) = 1.0 / odzq(k);
          ot(k) = pow(opres(k) * 1.e-5, rd_inv_cp) * oth(k);
          orho(k) = opres(k) / (rd * ot(k));
          oinv_rho(k) = 1.0 / orho(k);
          orhofacr(k) = pow(rhosur * oinv_rho(k), 0.54);
        });
      team.team_barrier();

      bool log_qxpresent;
      const Int k_qxtop = find_top(team, osqr, qsmall, kbot, ktop, kdir, log_qxpresent);

      if (log_qxpresent) {
        Real dt_left = dt;    // time remaining for sedi over full model (mp) time step
        Real prt_accum = 0.0; // precip rate for individual category

        Int k_qxbot = find_bottom(team, osqr, qsmall, kbot, k_qxtop, kdir, log_qxpresent);

        while (dt_left > 1.e-4) {
          Real Co_max = 0.0;
          Int kmin, kmax;

          Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, m.num_pack), [&] (Int k) {
              oV_qr(k) = 0;
              oV_nr(k) = 0;
            });
          team.team_barrier();

          util::set_min_max(k_qxbot, k_qxtop, kmin, kmax, RealSmallPack::n);
          Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, kmax-kmin+1), [&] (int pk_, Real& lmax) {
              const int pk = kmin + pk_;
              auto qr_gt_small = (olqr(pk) > qsmall);
              if (qr_gt_small.any()) {
                // Compute Vq, Vn:
                olnr(pk).set(qr_gt_small, max(olnr(pk), nsmall));
                Table3 t;
                RealSmallPack tmp1, tmp2;
                get_rain_dsd2_kokkos(qr_gt_small, olqr(pk), olnr(pk), olmu_r(pk),
                                     t.rdumii, t.dumii, ollamr(pk),
                                     m.mu_r_table, tmp1, tmp2);
                find_lookupTable_indices_3_kokkos(qr_gt_small, t, olmu_r(pk), ollamr(pk));
                // mass-weighted fall speed:
                olV_qr(pk).set(qr_gt_small,
                               apply_table(qr_gt_small, m.vm_table, t) * olrhofacr(pk));
                // number-weighted fall speed:
                olV_nr(pk).set(qr_gt_small,
                               apply_table(qr_gt_small, m.vn_table, t) * olrhofacr(pk));
                const auto Co_max_local = max(qr_gt_small, -1,
                                              olV_qr(pk) * dt_left * olinv_dzq(pk));
                if (Co_max_local > lmax)
                  lmax = Co_max_local;
              }
            }, Kokkos::Max<Real>(Co_max));
          team.team_barrier();
          if (Co_max < 0) {
            // qr is everywhere too small. Exit dt_left loop.
            break;
          }

          // compute dt_sub
          const Int Co_max_p1 = static_cast<Int>(Co_max + 1.0);
          const Real dt_sub = util::min(dt_left, dt_left / Co_max_p1);

          // Move bottom cell down by 1 if not at ground already.
          const Int k_temp = (k_qxbot == kbot) ? k_qxbot : k_qxbot - kdir;

          calc_first_order_upwind_step<2>(
            olrho, olinv_rho, olinv_dzq, team,
            m.num_vert, k_temp, k_qxtop, kdir, dt_sub,
            {&olflux_qx, &olflux_nx}, {&olV_qr, &olV_nr}, {&olqr, &olnr});
          team.team_barrier();

          // accumulated precip during time step
          if (k_qxbot == kbot) prt_accum += osflux_qx(kbot) * dt_sub;

          dt_left -= dt_sub;  // update time remaining for sedimentation
          if (k_qxbot != kbot) k_qxbot -= kdir;
        }

        Kokkos::single(
          Kokkos::PerTeam(team), [&] () {
            prt_liq(i) += prt_accum * Globals<Real>::INV_RHOW * odt;
          });
      }
    });
}

} // namespace p3
} // namespace micro_sed

#endif
