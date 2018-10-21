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

  using pack_t = RealPack;

  static constexpr const char* NAME = "kokkos_packnoiws";

  int num_horz, num_vert, num_pack;

  kokkos_2d_table_t<Real> vn_table, vm_table;
  kokkos_1d_table_t<Real> mu_r_table;

  team_policy policy;
  util::WorkspaceManager<RealPack> workspace_mgr;

public:
  MicroSedFuncPackNoiWsKokkos(int num_horz_, int num_vert_) :
    num_horz(num_horz_), num_vert(num_vert_),
    num_pack(scream::pack::npack<RealPack>(num_vert_)),
    vn_table("VN_TABLE"), vm_table("VM_TABLE"),
    mu_r_table("MU_R_TABLE"),
    policy(util::ExeSpaceUtils<>::get_default_team_policy(num_horz, num_pack)),
    workspace_mgr(num_pack, 11, policy) // rain sed's high-water is 11 workspace for any team
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

void micro_sed_func (
  MicroSedFuncPackNoiWsKokkos<Real>& m,
  const Int kts, const Int kte, const int its, const int ite, const Real dt,
  const kokkos_2d_t<RealPack>& qr, const kokkos_2d_t<RealPack>& nr,
  const kokkos_2d_t<RealPack>& th, const kokkos_2d_t<RealPack>& dzq, const kokkos_2d_t<RealPack>& pres,
  const kokkos_1d_t<Real>& prt_liq)
{
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

      auto workspace = m.workspace_mgr.get_workspace(team);

      Unmanaged<kokkos_1d_t<RealSmallPack> > linv_dzq, lrho, linv_rho, lrhofacr, lt;
      {
        Kokkos::Array<Unmanaged<kokkos_1d_t<RealSmallPack> >*, 5> ptrs = { {&linv_dzq, &lrho, &linv_rho, &lrhofacr, &lt} };
        Kokkos::Array<const char*, 5> names = { {"inv_dzq", "rho", "inv_rho", "rhofacr", "t"} };

        workspace.take_many(names, ptrs);
      }

      {
        const Unmanaged<kokkos_1d_t<RealPack> >
          odzq     = util::subview(dzq, i),
          oth      = util::subview(th, i),
          opres    = util::subview(pres, i);

        auto inv_dzq = biggize(linv_dzq);
        auto rho     = biggize(lrho);
        auto inv_rho = biggize(linv_rho);
        auto rhofacr = biggize(lrhofacr);
        auto t       = biggize(lt);

        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, m.num_pack), [&] (Int k) {
            inv_dzq(k) = 1.0 / odzq(k);
            t(k) = pow(opres(k) * 1.e-5, rd_inv_cp) * oth(k);
            rho(k) = opres(k) / (rd * t(k));
            inv_rho(k) = 1.0 / rho(k);
            rhofacr(k) = pow(rhosur * inv_rho(k), 0.54);
          });
        team.team_barrier();
      }

      const Unmanaged<kokkos_1d_t<Real> >
        osqr = util::subview(scalarize(qr), i);

      bool log_qxpresent;
      const Int k_qxtop = find_top(team, osqr, qsmall, kbot, ktop, kdir, log_qxpresent);

      if (log_qxpresent) {
        Real dt_left = dt;    // time remaining for sedi over full model (mp) time step
        Real prt_accum = 0.0; // precip rate for individual category

        Int k_qxbot = find_bottom(team, osqr, qsmall, kbot, k_qxtop, kdir, log_qxpresent);

        Unmanaged<kokkos_1d_t<RealSmallPack> > lV_qr, lV_nr, lflux_qx, lflux_nx, lmu_r, llamr;
        {
          Kokkos::Array<Unmanaged<kokkos_1d_t<RealSmallPack> >*, 6> ptrs = { {&lV_qr, &lV_nr, &lflux_qx, &lflux_nx, &lmu_r, &llamr} };
          Kokkos::Array<const char*, 6> names = { {"V_qr", "V_nr", "flux_qx", "flux_nx", "mu_r", "lamr"} };

          workspace.take_many(names, ptrs);
        }

        while (dt_left > 1.e-4) {
          Real Co_max = 0.0;
          Int kmin, kmax;

          const Unmanaged<kokkos_1d_t<RealSmallPack> >
            olqr = util::subview(smallize(qr), i),
            olnr = util::subview(smallize(nr), i);

          {
            auto V_qr = biggize(lV_qr);
            auto V_nr = biggize(lV_nr);
            Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, m.num_pack), [&] (Int k) {
                V_qr(k) = 0;
                V_nr(k) = 0;
              });
            team.team_barrier();
          }

          util::set_min_max(k_qxbot, k_qxtop, kmin, kmax, RealSmallPack::n);

          Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, kmax-kmin+1), [&] (int pk_, Real& lmax) {
              const int pk = kmin + pk_;
              auto qr_gt_small = (olqr(pk) > qsmall);
              if (qr_gt_small.any()) {
                // Compute Vq, Vn:
                olnr(pk).set(qr_gt_small, max(olnr(pk), nsmall));
                Table3 table;
                RealSmallPack tmp1, tmp2;
                get_rain_dsd2_kokkos(qr_gt_small, olqr(pk), olnr(pk), lmu_r(pk),
                                     table.rdumii, table.dumii, llamr(pk),
                                     m.mu_r_table, tmp1, tmp2);
                find_lookupTable_indices_3_kokkos(qr_gt_small, table, lmu_r(pk), llamr(pk));
                // mass-weighted fall speed:
                lV_qr(pk).set(qr_gt_small,
                              apply_table(qr_gt_small, m.vm_table, table) * lrhofacr(pk));
                // number-weighted fall speed:
                lV_nr(pk).set(qr_gt_small,
                              apply_table(qr_gt_small, m.vn_table, table) * lrhofacr(pk));
                const auto Co_max_local = max(qr_gt_small, -1,
                                              lV_qr(pk) * dt_left * linv_dzq(pk));
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
            lrho, linv_rho, linv_dzq, team,
            m.num_vert, k_temp, k_qxtop, kdir, dt_sub,
            {&lflux_qx, &lflux_nx}, {&lV_qr, &lV_nr}, {&olqr, &olnr});
          team.team_barrier();

          // accumulated precip during time step
          {
            auto sflux_qx = scalarize(lflux_qx);
            if (k_qxbot == kbot) prt_accum += sflux_qx(kbot) * dt_sub;
          }

          dt_left -= dt_sub;  // update time remaining for sedimentation
          if (k_qxbot != kbot) k_qxbot -= kdir;
        }

        Kokkos::single(
          Kokkos::PerTeam(team), [&] () {
            prt_liq(i) += prt_accum * Globals<Real>::INV_RHOW * odt;
          });
      }
      // Since we're at the end of the kernel, call reset instead of
      // individual releases.
      workspace.reset();
    });

  // m.workspace_mgr.report(); // uncomment for detailed debug info
}

} // namespace p3
} // namespace micro_sed

#endif
