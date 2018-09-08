#ifndef MICRO_SED_PACK_KOKKOS_HPP
#define MICRO_SED_PACK_KOKKOS_HPP

#include "util.hpp"
#include "initial_conditions.hpp"
#include "micro_kokkos.hpp"
#include "micro_sed_vanilla.hpp"
#include "scream_pack.hpp"

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>

namespace p3 {
namespace micro_sed_pack {

using micro_sed_vanilla::Globals;
using namespace scream;

#define SCREAM_PACKN 16
using RealPack = Pack<Real, SCREAM_PACKN>;
using IntPack = Pack<Int, SCREAM_PACKN>;

KOKKOS_FORCEINLINE_FUNCTION
kokkos_2d_t<Real> scalarize (const kokkos_2d_t<RealPack>& vp) {
  return Kokkos::View<Real**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >(
    reinterpret_cast<Real*>(vp.data()), vp.extent_int(0), SCREAM_PACKN * vp.extent_int(1));
}

// NOT for general use. This is just for dev work.
KOKKOS_FORCEINLINE_FUNCTION
kokkos_2d_t<RealPack> packize (const kokkos_2d_t<Real>& vp) {
  assert(vp.extent_int(1) % SCREAM_PACKN == 0);
  return Kokkos::View<RealPack**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >(
    reinterpret_cast<RealPack*>(vp.data()), vp.extent_int(0), vp.extent_int(1) / SCREAM_PACKN);
}

template <typename Real>
using kokkos_2d_table_t = Kokkos::View<Real[300][10], Layout, MemSpace>;

template <typename Real>
using kokkos_1d_table_t = Kokkos::View<Real[150], Layout, MemSpace>;

struct Table3 {
  int dumii, dumjj;
  Real rdumii, rdumjj, inv_dum3;
};

// Finds indices in rain lookup table (3)
KOKKOS_INLINE_FUNCTION
void find_lookupTable_indices_3_kokkos (Table3& t, const Real& mu_r, const Real& lamr)
{
  // find location in scaled mean size space
  const auto dum1 = (mu_r+1.) / lamr;
  if (dum1 <= 195.e-6) {
    t.inv_dum3  = 0.1;
    t.rdumii = (dum1*1.e6+5.)*t.inv_dum3;
    t.rdumii = util::max<Real>(t.rdumii, 1.);
    t.rdumii = util::min<Real>(t.rdumii,20.);
    t.dumii  = static_cast<int>(t.rdumii);
    t.dumii  = util::max(t.dumii, 1);
    t.dumii  = util::min(t.dumii,20);
  }
  else {
    t.inv_dum3  = Globals<Real>::THRD*0.1;
    t.rdumii = (dum1*1.e+6-195.)*t.inv_dum3 + 20.;
    t.rdumii = util::max<Real>(t.rdumii, 20.);
    t.rdumii = util::min<Real>(t.rdumii,300.);
    t.dumii  = static_cast<int>(t.rdumii);
    t.dumii  = util::max(t.dumii, 20);
    t.dumii  = util::min(t.dumii,299);
  }

  // find location in mu_r space
  t.rdumjj = mu_r+1.;
  t.rdumjj = util::max<Real>(t.rdumjj,1.);
  t.rdumjj = util::min<Real>(t.rdumjj,10.);
  t.dumjj  = static_cast<int>(t.rdumjj);
  t.dumjj  = util::max(t.dumjj,1);
  t.dumjj  = util::min(t.dumjj,9);
}

KOKKOS_INLINE_FUNCTION
Real apply_table (const kokkos_2d_table_t<Real>& table, const Table3& t) {
  const auto dum1 = (table(t.dumii-1, t.dumjj-1) + (t.rdumii-t.dumii) * t.inv_dum3 *
                     (table(t.dumii, t.dumjj-1) - table(t.dumii-1, t.dumjj-1)));
  const auto dum2 = (table(t.dumii-1, t.dumjj) + (t.rdumii-t.dumii) * t.inv_dum3 *
                     (table(t.dumii, t.dumjj) - table(t.dumii-1, t.dumjj)));
  return dum1 + (t.rdumjj - t.dumjj) * (dum2 - dum1);
}

// Computes and returns rain size distribution parameters
KOKKOS_INLINE_FUNCTION
void get_rain_dsd2_kokkos (
  const Real& qr, Real& nr, Real& mu_r, Real& rdumii, int& dumii, Real& lamr,
  const kokkos_1d_table_t<Real>& mu_r_table, Real& cdistr, Real& logn0r)
{
  constexpr Real nsmall = Globals<Real>::NSMALL;
  if (qr >= Globals<Real>::QSMALL) {
    // use lookup table to get mu
    // mu-lambda relationship is from Cao et al. (2008), eq. (7)

    // find spot in lookup table
    // (scaled N/q for lookup table parameter space_
    nr = util::max(nr, nsmall);
    const auto inv_dum = std::pow(qr / (Globals<Real>::CONS1 * nr * 6.0), Globals<Real>::THRD);

    if (inv_dum < 282.e-6) {
      mu_r = 8.282;
    }
    else if (inv_dum >= 282.e-6 && inv_dum < 502.e-6) {
      // interpolate
      rdumii = (inv_dum-250.e-6)*1.e+6*0.5;
      rdumii = util::max<Real>(rdumii,1.0);
      rdumii = util::min<Real>(rdumii,150.0);
      dumii  = static_cast<int>(rdumii);
      dumii  = util::min(149,dumii);
      mu_r   = mu_r_table(dumii-1) + (mu_r_table(dumii) - mu_r_table(dumii-1)) * (rdumii-dumii);
    }
    else if (inv_dum >= 502.e-6) {
      mu_r = 0.0;
    }

    lamr   = std::pow((Globals<Real>::CONS1 *nr *(mu_r+3.0) * (mu_r+2) * (mu_r+1.)/(qr)),
                      Globals<Real>::THRD); // recalculate slope based on mu_r
    const auto lammax = (mu_r+1.)*1.e+5;  // check for slope
    const auto lammin = (mu_r+1.)*1250.0; // set to small value since breakup is explicitly included (mean size 0.8 mm)

    // apply lambda limiters for rain
    const auto lt = lamr < lammin;
    const auto gt = lamr > lammax;
    if (lt || gt) {
      lamr = lt ? lammin : lammax;
      nr   = std::exp(3.*std::log(lamr) + std::log(qr) +
                      std::log(std::tgamma(mu_r+1.)) - std::log(std::tgamma(mu_r+4.)))/
        (Globals<Real>::CONS1);
    }

    cdistr  = nr/std::tgamma(mu_r+1.);
    // note: logn0r is calculated as log10(n0r);
    logn0r  = std::log10(nr) + (mu_r+1.)*std::log10(lamr) - std::log10(std::tgamma(mu_r+1));
  }
  else {
    lamr   = 0.0;
    cdistr = 0.0;
    logn0r = 0.0;
  }
}

struct MicroSedFuncPackKokkos {
  int num_horz, num_vert;

  // re-usable scratch views
  kokkos_2d_t<Real> V_qr, V_nr, flux_qx, flux_nx, mu_r, lamr, rhofacr, inv_dzq, rho, inv_rho, t, tmparr1;
  kokkos_2d_table_t<Real> vn_table, vm_table;
  kokkos_1d_table_t<Real> mu_r_table;

public:
  MicroSedFuncPackKokkos(int num_horz_, int num_vert_) :
    num_horz(num_horz_), num_vert(num_vert_),
    V_qr("V_qr", num_horz, num_vert),
    V_nr("V_nr", num_horz, num_vert),
    flux_qx("flux_qx", num_horz, num_vert),
    flux_nx("flux_nx", num_horz, num_vert),
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

void reset (MicroSedFuncPackKokkos& msvk) {
  Kokkos::parallel_for(
    "2d reset",
    util::ExeSpaceUtils<>::get_default_team_policy(msvk.num_horz, msvk.num_vert),
    KOKKOS_LAMBDA(const member_type& team_member) {
      const int i = team_member.league_rank();
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team_member, msvk.num_vert), [&] (int k) {
          msvk.V_qr(i, k)    = 0.0;
          msvk.V_nr(i, k)    = 0.0;
          msvk.flux_qx(i, k) = 0.0;
          msvk.flux_nx(i, k) = 0.0;
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

KOKKOS_INLINE_FUNCTION
void calc_first_order_upwind_step (
  const MicroSedFuncPackKokkos& m, const member_type& team_member, const int& i,
  const int& k_temp, const int& k_qxtop, const int& kdir, const Real& dt_sub,
  const kokkos_2d_t<Real>& flux, const kokkos_2d_t<Real>& V, const kokkos_2d_t<Real>& r)
{
  int kmin, kmax;

  // calculate fluxes
  util::set_min_max(k_temp, k_qxtop + kdir, kmin, kmax);
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team_member, kmax-kmin+1), [&] (int k_) {
      const int k = kmin + k_;
      flux(i, k) = V(i, k) * r(i, k) * m.rho(i, k);
    });
  team_member.team_barrier();

  Kokkos::single(
    Kokkos::PerTeam(team_member), [&] () {
      // for top level only (since flux is 0 above)
      int k = k_qxtop;
      // compute flux divergence
      const auto fluxdiv = -flux(i, k) * m.inv_dzq(i, k);
      // update prognostic variables
      r(i, k) += fluxdiv * dt_sub * m.inv_rho(i, k);
    });

  util::set_min_max(k_qxtop - kdir, k_temp, kmin, kmax);
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team_member, kmax-kmin+1), [&] (int k_) {
      const int k = kmin + k_;
      // compute flux divergence
      const auto fluxdiv = (flux(i, k+kdir) - flux(i, k)) * m.inv_dzq(i, k);
      // update prognostic variables
      r(i, k) += fluxdiv * dt_sub * m.inv_rho(i, k);
    });
}

void micro_sed_func_pack_kokkos (
  MicroSedFuncPackKokkos& m,
  const int kts, const int kte, const int its, const int ite, const Real dt,
  const kokkos_2d_t<Real>& qr, const kokkos_2d_t<Real>& nr,
  const kokkos_2d_t<Real>& th, const kokkos_2d_t<Real>& dzq, const kokkos_2d_t<Real>& pres,
  const kokkos_1d_t<Real>& prt_liq)
{
  // constants
  const Real odt = 1.0 / dt;
  constexpr Real nsmall = Globals<Real>::NSMALL;

  // direction of vertical leveling
  const int kbot = (kts < kte) ? 0: m.num_vert-1;
  const int kdir = (kts < kte) ? 1  : -1;

  const kokkos_2d_t<RealPack>
    pdzq = packize(dzq),
    ppres = packize(pres),
    pinv_dzq = packize(m.inv_dzq),
    pt = packize(m.t),
    pth = packize(th),
    prho = packize(m.rho),
    pinv_rho = packize(m.inv_rho),
    prhofacr = packize(m.rhofacr);

  // Rain sedimentation:  (adaptivive substepping)
  Kokkos::parallel_for(
    "main rain sed loop",
    util::ExeSpaceUtils<>::get_default_team_policy(m.num_horz, m.num_vert),
    KOKKOS_LAMBDA(const member_type& team_member) {
      const int i = team_member.league_rank();

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team_member, m.num_vert / SCREAM_PACKN), [&] (int k) {
          // inverse of thickness of layers
          pinv_dzq(i, k) = 1 / pdzq(i, k);
          pt(i, k) = pow(ppres(i, k) * 1.e-5, Globals<Real>::RD * Globals<Real>::INV_CP) * pth(i, k);
          prho(i, k) = ppres(i, k) / (Globals<Real>::RD * pt(i, k));
          pinv_rho(i, k) = 1.0 / prho(i, k);
          prhofacr(i, k) = pow(Globals<Real>::RHOSUR * pinv_rho(i, k), 0.54);
        });
      team_member.team_barrier();

      // Note, we are skipping supersaturation checks

      bool log_qxpresent = false;
      int k_qxtop = -1; // avoid warning, but don't use a meaningful value

      // find top, determine qxpresent
      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team_member, m.num_vert), [&] (int k, int& lmax) {
          if (qr(i, k) >= Globals<Real>::QSMALL && k*kdir > lmax) {
            lmax = k*kdir;
          }
        }, Kokkos::Max<int>(k_qxtop));
      // If the if statement in the parallel_reduce is never true, k_qxtop will
      // end up being a large negative number, Max::init()'s value.
      k_qxtop *= kdir;
      log_qxpresent = k_qxtop >= 0;

      if (log_qxpresent) {
        Real dt_left = dt;    // time remaining for sedi over full model (mp) time step
        Real prt_accum = 0.0; // precip rate for individual category
        int k_qxbot = 0;

        // find bottom
        Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team_member, m.num_vert), [&] (int k, int& lmin) {
            if (qr(i, k) >= Globals<Real>::QSMALL && k*kdir < lmin) {
              lmin = k*kdir;
            }
          }, Kokkos::Min<int>(k_qxbot));
        // As log_qxpresent is true, we don't have to worry about this reduction
        // as we did for the one for k_qxtop.
        k_qxbot *= kdir;

        while (dt_left > 1.e-4) {
          Real Co_max = 0.0;
          Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_member, m.num_vert),
            [&] (int kk) {
              m.V_qr(i, kk) = 0.0;
              m.V_nr(i, kk) = 0.0;
            });
          team_member.team_barrier();

          int kmin, kmax;
          util::set_min_max(k_qxtop, k_qxbot, kmin, kmax);
          Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team_member, kmax-kmin+1), [&] (int k_, Real& lmax) {
              const int k = kmin + k_;
              if (qr(i, k) > Globals<Real>::QSMALL) {
                // Compute Vq, Vn:
                nr(i, k) = util::max(nr(i, k), nsmall);
                Real tmp1, tmp2;
                Table3 t;
                get_rain_dsd2_kokkos(qr(i, k), nr(i, k), m.mu_r(i, k), t.rdumii, t.dumii,
                                     m.lamr(i, k), m.mu_r_table, tmp1, tmp2);
                find_lookupTable_indices_3_kokkos(t, m.mu_r(i, k), m.lamr(i, k));
                // mass-weighted fall speed:
                m.V_qr(i, k) = apply_table(m.vm_table, t) * m.rhofacr(i, k);
                // number-weighted fall speed:
                m.V_nr(i, k) = apply_table(m.vn_table, t) * m.rhofacr(i, k);
              }
              Real Co_max_local = m.V_qr(i, k) * dt_left * m.inv_dzq(i, k);
              if (Co_max_local > lmax)
                lmax = Co_max_local;
            }, Kokkos::Max<Real>(Co_max));

          // compute dt_sub
          const int tmpint1 = static_cast<int>(Co_max + 1.0);
          const Real dt_sub = util::min(dt_left, dt_left / tmpint1);

          const int k_temp = (k_qxbot == kbot) ? k_qxbot : (k_qxbot - kdir);

          calc_first_order_upwind_step(m, team_member, i,
                                       k_temp, k_qxtop, kdir, dt_sub,
                                       m.flux_qx, m.V_qr, qr);
          calc_first_order_upwind_step(m, team_member, i,
                                       k_temp, k_qxtop, kdir, dt_sub,
                                       m.flux_nx, m.V_nr, nr);


          // accumulated precip during time step
          if (k_qxbot == kbot)
            prt_accum += m.flux_qx(i, kbot) * dt_sub;

          dt_left -= dt_sub;  // update time remaining for sedimentation
          if (k_qxbot != kbot)
            k_qxbot -= kdir;
        }

        Kokkos::single(
          Kokkos::PerTeam(team_member), [&]() {
            prt_liq(i) += prt_accum * Globals<Real>::INV_RHOW * odt;
          });
      }
    });

  reset(m);
}

void populate_kokkos_from_vec (
  const int num_horz, const int num_vert, vector_2d_t<Real> const& vec, kokkos_2d_t<Real>& device)
{
  typename kokkos_2d_t<Real>::HostMirror mirror = Kokkos::create_mirror_view(device);

  for (int i = 0; i < num_horz; ++i) {
    for (int k = 0; k < num_vert; ++k) {
      mirror(i, k) = vec[i][k];
    }
  }

  Kokkos::deep_copy(device, mirror);
}

void dump_to_file_k (
  const kokkos_2d_t<Real>& qr, const kokkos_2d_t<Real>& nr, const kokkos_2d_t<Real>& th,
  const kokkos_2d_t<Real>& dzq, const kokkos_2d_t<Real>& pres, const kokkos_1d_t<Real>& prt_liq,
  const Real dt, const int ts)
{
  const int ni = qr.extent_int(0);
  const int nk = qr.extent_int(1);

  auto qr_m      = Kokkos::create_mirror_view(qr);
  auto nr_m      = Kokkos::create_mirror_view(nr);
  auto th_m      = Kokkos::create_mirror_view(th);
  auto dzq_m     = Kokkos::create_mirror_view(dzq);
  auto pres_m    = Kokkos::create_mirror_view(pres);
  auto prt_liq_m = Kokkos::create_mirror_view(prt_liq);

  Kokkos::deep_copy(qr_m,      qr);
  Kokkos::deep_copy(nr_m,      nr);
  Kokkos::deep_copy(th_m,      th);
  Kokkos::deep_copy(dzq_m,     dzq);
  Kokkos::deep_copy(pres_m,    pres);
  Kokkos::deep_copy(prt_liq_m, prt_liq);

  p3::micro_sed_vanilla::dump_to_file(
    "pack_kokkos", qr_m.data(), nr_m.data(), th_m.data(), dzq_m.data(), pres_m.data(),
    prt_liq_m.data(), ni, nk, dt, ts);
}

void micro_sed_func_pack_kokkos_wrap (
  const int ni, const int nk, const Real dt, const int ts, const int kdir)
{
  std::vector<Real> v(nk);
  vector_2d_t<Real> qr_v(ni, v), nr_v(ni, v), th_v(ni, v), dzq_v(ni, v), pres_v(ni, v);

  util::dump_arch();
  std::cout << "Running micro_sed_pack_kokkos with ni=" << ni << ", nk=" << nk
            << ", dt=" << dt << ", ts=" << ts << ", kdir=" << kdir << std::endl;

  p3::micro_sed_vanilla::populate_input(ni, nk, kdir, qr_v, nr_v, th_v, dzq_v, pres_v);

  MicroSedFuncPackKokkos mspk(ni, nk);
  const int np = mspk.num_vert;
  
  kokkos_2d_t<Real> qr("qr", ni, np), nr("nr", ni, np), th("th", ni, np), dzq("dzq", ni, np),
    pres("pres", ni, np);
  kokkos_1d_t<Real> prt_liq("prt_liq", ni);

  for (auto item : { std::make_pair(&qr_v, &qr), std::make_pair(&nr_v, &nr), std::make_pair(&th_v, &th),
        std::make_pair(&dzq_v, &dzq), std::make_pair(&pres_v, &pres)}) {
    populate_kokkos_from_vec(ni, nk, *(item.first), *(item.second));
  }

  auto start = std::chrono::steady_clock::now();

  for (int i = 0; i < ts; ++i) {
    micro_sed_func_pack_kokkos(mspk,
                               kdir == 1 ? 1 : nk, kdir == 1 ? nk : 1,
                               1, ni, dt, qr, nr, th, dzq, pres, prt_liq);
    Kokkos::fence();
  }

  auto finish = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
  printf("Time = %1.3e seconds\n", 1e-6*duration.count());

  dump_to_file_k(qr, nr, th, dzq, pres, prt_liq, dt, ts);
}

} // namespace p3
} // namespace micro_sed_pack

#endif
