#ifndef MICRO_SED_WORKSPACE_KOKKOS_HPP
#define MICRO_SED_WORKSPACE_KOKKOS_HPP

#include "util.hpp"
#include "initial_conditions.hpp"
#include "micro_kokkos.hpp"
#include "micro_sed_vanilla_kokkos.hpp"

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>

namespace p3 {
namespace micro_sed {

template <typename Real>
struct MicroSedFuncWorkspaceKokkos
{
  int num_horz, num_vert, concurrency;

  //
  // re-usable scratch views
  //
  kokkos_2d_t<Real> V_qr, V_nr, flux_qx, flux_nx; // 2d to prevent race conditions
  kokkos_2d_t<Real> mu_r, lamr, rhofacr, inv_dzq, rho, inv_rho, t, tmparr1; // legit 2d

  kokkos_2d_table_t<Real> vn_table, vm_table;
  kokkos_1d_table_t<Real> mu_r_table;

  static constexpr char* NAME = "kokkos_workspace";

public:
  MicroSedFuncWorkspaceKokkos(int num_horz_, int num_vert_) :
    num_horz(num_horz_), num_vert(num_vert_),
    concurrency(util::ExeSpaceUtils<>::get_num_concurrent_teams(util::ExeSpaceUtils<>::get_default_team_policy(num_horz, num_vert))),
    V_qr("V_qr", concurrency, num_vert),
    V_nr("V_nr", concurrency, num_vert),
    flux_qx("flux_qx", concurrency, num_vert),
    flux_nx("flux_nx", concurrency, num_vert),
    mu_r("mu_r", num_horz, num_vert),
    lamr("lamr", num_horz, num_vert),
    rhofacr("rhofacr", num_horz, num_vert),
    inv_dzq("inv_dzq", num_horz, num_vert),
    rho("rho", num_horz, num_vert),
    inv_rho("inv_rho", num_horz, num_vert),
    t("t", num_horz, num_vert),
    tmparr1("tmparr1", num_horz, num_vert),
    vn_table("VN_TABLE"), vm_table("VM_TABLE"),
    mu_r_table("MU_R_TABLE")
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
};

template <typename Real>
void reset(MicroSedFuncWorkspaceKokkos<Real>& msvk)
{
  Kokkos::parallel_for("2d reset",
                       util::ExeSpaceUtils<>::get_default_team_policy(msvk.num_horz, msvk.num_vert),
                       KOKKOS_LAMBDA(member_type team_member) {
    const int i = team_member.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int k) {
      if (i < msvk.concurrency) {
        msvk.V_qr(i, k)    = 0.0;
        msvk.V_nr(i, k)    = 0.0;
        msvk.flux_qx(i, k) = 0.0;
        msvk.flux_nx(i, k) = 0.0;
      }
      msvk.mu_r(i, k)    = 0.0;
      msvk.lamr(i, k)    = 0.0;
      msvk.rhofacr(i, k) = 0.0;
      msvk.inv_dzq(i, k) = 0.0;
      msvk.rho(i, k)     = 0.0;
      msvk.inv_rho(i, k) = 0.0;
      msvk.t(i, k)       = 0.0;
      msvk.tmparr1(i, k) = 0.0;
    });
  });
}

/**
 * Arg explanation
 *
 * kts: vertical array bound (top)
 * kte: vertical array bound (bottom)
 * ni: number of columns in slab
 * nk: number of vertical levels
 * its: horizontal array bound
 * ite: horizontal array bound
 * dt: time step
 * qr: rain, mass mixing ratio  (in/out)
 * nr: rain, number mixing ratio (in/out)
 * th: potential temperature                    K
 * dzq: vertical grid spacing                   m
 * pres: pressure                               Pa
 * prt_liq: precipitation rate, total liquid    m s-1  (output)
 */
void micro_sed_func(MicroSedFuncWorkspaceKokkos<Real>& msvk,
                    const int kts, const int kte, const int its, const int ite, const Real dt,
                    kokkos_2d_t<Real> & qr, kokkos_2d_t<Real> & nr,
                    kokkos_2d_t<Real> const& th, kokkos_2d_t<Real> const& dzq, kokkos_2d_t<Real> const& pres,
                    kokkos_1d_t<Real> & prt_liq)
{
  // constants
  const Real odt = 1.0 / dt;
  constexpr Real nsmall = Globals<Real>::NSMALL;

  // direction of vertical leveling
#ifdef TRACE
  const int ktop = (kts < kte) ? msvk.num_vert-1 : 0;
#endif
  const int kbot = (kts < kte) ? 0: msvk.num_vert-1;
  const int kdir = (kts < kte) ? 1  : -1;

  // Rain sedimentation:  (adaptivive substepping)
  trace_loop("i_loop_main", 0, msvk.num_horz);
  Kokkos::parallel_for("main rain sed loop",
                       util::ExeSpaceUtils<>::get_default_team_policy(msvk.num_horz, msvk.num_vert),
                       KOKKOS_LAMBDA(member_type team_member) {
    const int i = team_member.league_rank();
    const int i_team = i % msvk.concurrency;

    trace_loop("  k_loop_1", kbot, ktop);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int k) {
      // inverse of thickness of layers
      msvk.inv_dzq(i, k) = 1 / dzq(i, k);
      msvk.t(i, k) = std::pow(pres(i, k) * 1.e-5, Globals<Real>::RD * Globals<Real>::INV_CP) * th(i, k);
      msvk.rho(i, k) = pres(i, k) / (Globals<Real>::RD * msvk.t(i, k));
      msvk.inv_rho(i, k) = 1.0 / msvk.rho(i, k);
      msvk.rhofacr(i, k) = std::pow(Globals<Real>::RHOSUR * msvk.inv_rho(i, k), 0.54);
      trace_data("    rhofacr", i, k, msvk.rhofacr(i, k));
    });
    team_member.team_barrier();

    // Note, we are skipping supersaturation checks

    bool log_qxpresent = false;
    int k_qxtop = -1; // avoid warning, but don't use a meanigful value

    // find top, determine qxpresent
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int k, int& lmax) {
      if (qr(i, k) >= Globals<Real>::QSMALL && k*kdir > lmax) {
        lmax = k*kdir;
      }
    }, Kokkos::Max<int>(k_qxtop));
    // If the if statement in the parallel_reduce is never true,
    // k_qxtop will end up being a large negative number,
    // Max::init()'s value.
    k_qxtop *= kdir;
    log_qxpresent = k_qxtop >= 0;

    // JGF: It appears rain sedimentation is mostly nothing unless log_qxpresent is true
    if (log_qxpresent) {

      Real dt_left = dt;    // time remaining for sedi over full model (mp) time step
      Real prt_accum = 0.0; // precip rate for individual category
      int k_qxbot = 0;

      // find bottom
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int k, int& lmin) {
        if (qr(i, k) >= Globals<Real>::QSMALL && k*kdir < lmin) {
          lmin = k*kdir;
        }
      }, Kokkos::Min<int>(k_qxbot));
      // As log_qxpresent is true, we don't have to worry about this
      // reduction as we did for the one for k_qxtop.
      k_qxbot *= kdir;

      while (dt_left > 1.e-4) {
        Real Co_max = 0.0;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int kk) {
          msvk.V_qr(i_team, kk) = 0.0;
          msvk.V_nr(i_team, kk) = 0.0;
        });
        team_member.team_barrier();

        trace_loop("  k_loop_sedi_r1", k_qxtop, k_qxbot);
        int kmin, kmax;
        util::set_min_max(k_qxtop, k_qxbot, kmin, kmax);
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, kmax-kmin+1), [=] (int k_, Real& lmax) {
          const int k = kmin + k_;
          if (qr(i, k) > Globals<Real>::QSMALL) {
            // Compute Vq, Vn:
            nr(i, k) = util::max(nr(i, k), nsmall);
            trace_data("    nr", i, k, nr(i, k));
            Real rdumii=0.0, tmp1=0.0, tmp2=0.0, rdumjj=0.0, inv_dum3=0.0;
            int dumii=0, dumjj=0;
            get_rain_dsd2_kokkos(qr(i, k), nr(i, k), msvk.mu_r(i, k), rdumii, dumii, msvk.lamr(i, k), msvk.mu_r_table, tmp1, tmp2);
            find_lookupTable_indices_3_kokkos(dumii, dumjj, rdumii, rdumjj, inv_dum3, msvk.mu_r(i, k), msvk.lamr(i, k));

            // mass-weighted fall speed:
            Real dum1 = msvk.vm_table(dumii-1, dumjj-1) + (rdumii-dumii) * inv_dum3 *
              (msvk.vm_table(dumii, dumjj-1) - msvk.vm_table(dumii-1, dumjj-1));
            Real dum2 = msvk.vm_table(dumii-1, dumjj) + (rdumii-dumii) * inv_dum3 *
              (msvk.vm_table(dumii, dumjj) - msvk.vm_table(dumii-1, dumjj));

            msvk.V_qr(i_team, k) = (dum1 + (rdumjj - dumjj) * (dum2 - dum1)) * msvk.rhofacr(i, k);
            trace_data("    V_qr", i_team, k, msvk.V_qr(i_team, k));

            // number-weighted fall speed:
            dum1 = msvk.vn_table(dumii-1, dumjj-1) + (rdumii-dumii) * inv_dum3 *
              (msvk.vn_table(dumii, dumjj-1) - msvk.vn_table(dumii-1, dumjj-1));
            dum2 = msvk.vn_table(dumii-1, dumjj) + (rdumii-dumii) * inv_dum3 *
              (msvk.vn_table(dumii, dumjj) - msvk.vn_table(dumii-1, dumjj));

            msvk.V_nr(i_team, k) = (dum1 + (rdumjj - dumjj) * (dum2 - dum1)) * msvk.rhofacr(i, k);
            trace_data("    V_nr", i_team, k, msvk.V_nr(i_team, k));
          }
          Real Co_max_local = msvk.V_qr(i_team, k) * dt_left * msvk.inv_dzq(i, k);
          if (Co_max_local > lmax) {
            lmax = Co_max_local;
          }
          trace_data("  Co_max", 0, 0, Co_max);
        }, Kokkos::Max<Real>(Co_max));

        // compute dt_sub
        int tmpint1 = static_cast<int>(Co_max + 1.0);
        Real dt_sub = util::min(dt_left, dt_left / tmpint1);

        int k_temp = (k_qxbot == kbot) ? k_qxbot : (k_qxbot - kdir);

        // calculate fluxes
        trace_loop("  k_flux_loop", k_temp, k_qxtop);
        util::set_min_max(k_temp, k_qxtop+kdir, kmin, kmax);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, kmax-kmin+1), [&] (int k_) {
          const int k = kmin + k_;
          msvk.flux_qx(i_team, k) = msvk.V_qr(i_team, k) * qr(i, k) * msvk.rho(i, k);
          trace_data("    flux_qx", i_team, k, msvk.flux_qx(i_team, k));
          msvk.flux_nx(i_team, k) = msvk.V_nr(i_team, k) * nr(i, k) * msvk.rho(i, k);
          trace_data("    flux_nx", i_team, k, msvk.flux_nx(i_team, k));
        });
        team_member.team_barrier();

        // accumulated precip during time step
        if (k_qxbot == kbot) {
          prt_accum += msvk.flux_qx(i_team, kbot) * dt_sub;
        }

        Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
          // for top level only (since flux is 0 above)
          int k = k_qxtop;
          // compute flux divergence
          const Real fluxdiv_qx = -msvk.flux_qx(i_team, k) * msvk.inv_dzq(i, k);
          const Real fluxdiv_nx = -msvk.flux_nx(i_team, k) * msvk.inv_dzq(i, k);
          // update prognostic variables
          qr(i, k) += fluxdiv_qx * dt_sub * msvk.inv_rho(i, k);
          trace_data("  qr", i, k, qr(i, k));
          nr(i, k) += fluxdiv_nx * dt_sub * msvk.inv_rho(i, k);
          trace_data("  nr", i, k, nr(i, k));
        });
        team_member.team_barrier();

        trace_loop("  k_flux_div_loop", k_qxtop - kdir, k_temp);
        util::set_min_max(k_qxtop - kdir, k_temp, kmin, kmax);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, kmax-kmin+1), [&] (int k_) {
          const int k = kmin + k_;
          {
            // compute flux divergence
            const Real fluxdiv_qx = (msvk.flux_qx(i_team, k+kdir) - msvk.flux_qx(i_team, k)) * msvk.inv_dzq(i, k);
            const Real fluxdiv_nx = (msvk.flux_nx(i_team, k+kdir) - msvk.flux_nx(i_team, k)) * msvk.inv_dzq(i, k);
            // update prognostic variables
            qr(i, k) += fluxdiv_qx * dt_sub * msvk.inv_rho(i, k);
            trace_data("    qr", i, k, qr(i, k));
            nr(i, k) += fluxdiv_nx  *dt_sub * msvk.inv_rho(i, k);
            trace_data("    nr", i, k, nr(i, k));
          }
        });
        team_member.team_barrier();

        dt_left -= dt_sub;  // update time remaining for sedimentation
        if (k_qxbot != kbot) {
          k_qxbot -= kdir;
        }
      }

      trace_data("  prt_liq", i, 0, prt_liq(i));
      trace_data("  prt_accum", 0, 0, prt_accum);
      Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
        prt_liq(i) += prt_accum * Globals<Real>::INV_RHOW * odt;
      });
      trace_data("  prt_liq", i, 0, prt_liq(i));
    }
  });

  reset(msvk);
}

} // namespace p3
} // namespace micro_sed

#endif
