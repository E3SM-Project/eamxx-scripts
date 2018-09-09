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
#define SCREAM_SMALL_PACK_FACTOR 2
template <typename T> using BigPack = Pack<T, SCREAM_PACKN>;
template <typename T> using SmallPack = Pack<T, SCREAM_PACKN / SCREAM_SMALL_PACK_FACTOR>;
using RealPack = BigPack<Real>;
using IntPack = BigPack<Int>;
using RealSmallPack = SmallPack<Real>;
using IntSmallPack = SmallPack<Int>;

template <typename T> KOKKOS_FORCEINLINE_FUNCTION
kokkos_2d_t<T> scalarize (const kokkos_2d_t<BigPack<T> >& vp) {
  return Kokkos::View<T**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >(
    reinterpret_cast<T*>(vp.data()), vp.extent_int(0), RealPack::n * vp.extent_int(1));
}

// NOT for general use. This is just for dev work.
template <typename T> KOKKOS_FORCEINLINE_FUNCTION
kokkos_2d_t<BigPack<T> > packize (const kokkos_2d_t<T>& vp) {
  assert(vp.extent_int(1) % RealPack::n == 0);
  return Kokkos::View<BigPack<T>**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >(
    reinterpret_cast<BigPack<T>*>(vp.data()), vp.extent_int(0), vp.extent_int(1) / RealPack::n);
}
template <typename T> KOKKOS_FORCEINLINE_FUNCTION
kokkos_2d_t<SmallPack<T> > smallize (const kokkos_2d_t<T>& vp) {
  assert(vp.extent_int(1) % RealSmallPack::n == 0);
  return Kokkos::View<SmallPack<T>**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >(
    reinterpret_cast<SmallPack<T>*>(vp.data()), vp.extent_int(0), vp.extent_int(1) / RealSmallPack::n);
}

template <typename T> KOKKOS_FORCEINLINE_FUNCTION
kokkos_2d_t<SmallPack<T> > smallize (const kokkos_2d_t<BigPack<T> >& vp) {
  return Kokkos::View<SmallPack<T>**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >(
    reinterpret_cast<SmallPack<T>*>(vp.data()), vp.extent_int(0),
    SCREAM_SMALL_PACK_FACTOR * vp.extent_int(1));
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

    // These are not used in the micro-app, but keep them to model flops.
    cdistr = nr/std::tgamma(mu_r+1.);
    // note: logn0r is calculated as log10(n0r);
    logn0r = std::log10(nr) + (mu_r+1.)*std::log10(lamr) - std::log10(std::tgamma(mu_r+1));
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
  kokkos_2d_t<Real> V_qr, V_nr, flux, mu_r, lamr, rhofacr, inv_dzq, rho, inv_rho, t;
  kokkos_2d_table_t<Real> vn_table, vm_table;
  kokkos_1d_table_t<Real> mu_r_table;

public:
  MicroSedFuncPackKokkos(int num_horz_, int num_vert_) :
    num_horz(num_horz_), num_vert(num_vert_),
    V_qr("V_qr", num_horz, num_vert),
    V_nr("V_nr", num_horz, num_vert),
    flux("flux", num_horz, num_vert),
    mu_r("mu_r", num_horz, num_vert),
    lamr("lamr", num_horz, num_vert),
    rhofacr("rhofacr", num_horz, num_vert),
    inv_dzq("inv_dzq", num_horz, num_vert),
    rho("rho", num_horz, num_vert),
    inv_rho("inv_rho", num_horz, num_vert),
    t("t", num_horz, num_vert),
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

//TODO At the end, if it's still true that computing the steps for qr, nr in
// sequence is more expensive than together, take arrays of flux, V, r with
// compile-time size and iterate over them.

// Calculate the step in the region [k_bot, k_top].
template <int kdir>
KOKKOS_INLINE_FUNCTION
void calc_first_order_upwind_step (
  const MicroSedFuncPackKokkos& m, const member_type& team, const int& i,
  const int& k_bot, const int& k_top, const Real& dt_sub,
  const kokkos_2d_t<RealSmallPack>& flux, const kokkos_2d_t<RealSmallPack>& V,
  const kokkos_2d_t<RealSmallPack>& r)
{
  const kokkos_2d_t<RealSmallPack>
    rho = smallize(m.rho),
    inv_rho = smallize(m.inv_rho),
    inv_dzq = smallize(m.inv_dzq);

  int
    kmin = ( kdir == 1 ? k_bot : k_top)      / RealSmallPack::n,
    // Add 1 to make [kmin, kmax). But then the extra term (RealSmallPack::n -
    // 1) to determine pack index cancels the +1.
    kmax = ((kdir == 1 ? k_top : k_bot) + RealSmallPack::n) / RealSmallPack::n;

  // for top level only (since flux is 0 above)
  const int k_top_pack = k_top / RealSmallPack::n;
  flux(i, k_top_pack) = 0;
  team.team_barrier();

  // calculate fluxes
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, kmax - kmin), [&] (int k_) {
      const int k = kmin + k_;
      flux(i, k) = V(i, k) * r(i, k) * rho(i, k);
    });
  team.team_barrier();

  Kokkos::single(
    Kokkos::PerTeam(team), [&] () {
      const int k = k_top_pack;
      // compute flux divergence
      const auto flux_pkdir = (kdir == -1) ?
        shift_right(0, flux(i, k)) :
        shift_left (0, flux(i, k));
      const auto fluxdiv = (flux_pkdir - flux(i, k)) * inv_dzq(i, k);
      // update prognostic variables
      r(i, k) += fluxdiv * dt_sub * inv_rho(i, k);
    });

  if (kdir == 1)
    --kmax;
  else
    ++kmin;

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, kmax - kmin), [&] (int k_) {
      const int k = kmin + k_;
      // compute flux divergence
      const auto flux_pkdir = (kdir == -1) ?
        shift_right(flux(i, k+kdir), flux(i, k)) :
        shift_left (flux(i, k+kdir), flux(i, k));
      const auto fluxdiv = (flux_pkdir - flux(i, k)) * inv_dzq(i, k);
      // update prognostic variables
      r(i, k) += fluxdiv * dt_sub * inv_rho(i, k);
    });
}

KOKKOS_INLINE_FUNCTION
void calc_first_order_upwind_step (
  const MicroSedFuncPackKokkos& m, const member_type& team, const int& i,
  const int& k_bot, const int& k_top, const int& kdir, const Real& dt_sub,
  const kokkos_2d_t<RealSmallPack>& flux, const kokkos_2d_t<RealSmallPack>& V,
  const kokkos_2d_t<RealSmallPack>& r)
{
  if (kdir == 1)
    calc_first_order_upwind_step< 1>(
      m, team, i, k_bot, k_top, dt_sub, flux, V, r);
  else
    calc_first_order_upwind_step<-1>(
      m, team, i, k_bot, k_top, dt_sub, flux, V, r);
}

// Find the bottom and top of the mixing ratio, e.g., qr. It's worth casing
// these out in two ways: 1 thread/column vs many, and by kdir.
template <typename Array1>
KOKKOS_INLINE_FUNCTION
Int find_bottom (const member_type& team,
                 const Array1& v, const Real& small,
                 const Int& kbot, const Int& ktop, const Int& kdir,
                 bool& log_present) {
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
        Kokkos::TeamThreadRange(team, kbot - ktop + 1), [&] (int k_, int& lmax) {
          const int k = ktop + k_;
          if (v(k) >= small && k > lmax)
            lmax = k;
        }, Kokkos::Max<int>(k_xbot));
      log_present = k_xbot >= ktop;
    } else {
      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, ktop - kbot + 1), [&] (int k_, int& lmin) {
          const int k = kbot + k_;
          if (v(k) >= small && k < lmin)
            lmin = k;
        }, Kokkos::Min<int>(k_xbot));
      log_present = k_xbot <= ktop;
    }
  }
  return k_xbot;
}

template <typename Array1>
KOKKOS_INLINE_FUNCTION
Int find_top (const member_type& team,
              const Array1& v, const Real& small,
              const Int& kbot, const Int& ktop, const Int& kdir,
              bool& log_present) {
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
        Kokkos::TeamThreadRange(team, kbot - ktop + 1), [&] (int k_, int& lmin) {
          const int k = ktop + k_;
          if (v(k) >= small && k < lmin)
            lmin = k;
        }, Kokkos::Min<int>(k_xtop));
      log_present = k_xtop <= kbot;
    } else {
      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, ktop - kbot + 1), [&] (int k_, int& lmax) {
          const int k = kbot + k_;
          if (v(k) >= small && k > lmax)
            lmax = k;
        }, Kokkos::Max<int>(k_xtop));
      log_present = k_xtop >= kbot;
    }
  }
  return k_xtop;
}

template <typename Array1>
KOKKOS_INLINE_FUNCTION
Int find_bottom (const Array1& v, const Int& nk, const Int& kdir) {
  return -1;
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
  const int kbot = (kts < kte) ? 0 : m.num_vert-1;
  const int ktop = (kts < kte) ? m.num_vert-1 : 0;
  const int kdir = (kts < kte) ? 1 : -1;

  const kokkos_2d_t<RealPack>
    pdzq = packize(dzq),
    ppres = packize(pres),
    pinv_dzq = packize(m.inv_dzq),
    pt = packize(m.t),
    pth = packize(th),
    prho = packize(m.rho),
    pinv_rho = packize(m.inv_rho),
    prhofacr = packize(m.rhofacr),
    pV_qr = packize(m.V_qr),
    pV_nr = packize(m.V_nr),
    pflux = packize(m.flux),
    pqr = packize(qr),
    pnr = packize(nr);
  const kokkos_2d_t<RealSmallPack>
    sdzq = smallize(dzq),
    spres = smallize(pres),
    sinv_dzq = smallize(m.inv_dzq),
    st = smallize(m.t),
    sth = smallize(th),
    srho = smallize(m.rho),
    sinv_rho = smallize(m.inv_rho),
    srhofacr = smallize(m.rhofacr),
    sV_qr = smallize(m.V_qr),
    sV_nr = smallize(m.V_nr),
    sflux = smallize(m.flux),
    sqr = smallize(qr),
    snr = smallize(nr);

  const auto rd = Globals<Real>::RD;
  const auto rd_inv_cp = Globals<Real>::RD * Globals<Real>::INV_CP;
  const auto rhosur = Globals<Real>::RHOSUR;

  // Rain sedimentation:  (adaptive substepping)
  Kokkos::parallel_for(
    "main rain sed loop",
    util::ExeSpaceUtils<>::get_default_team_policy(m.num_horz, m.num_vert),
    KOKKOS_LAMBDA(const member_type& team) {
      const int i = team.league_rank();

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, m.num_vert / RealPack::n), [&] (int k) {
          // inverse of thickness of layers
          pinv_dzq(i, k) = 1 / pdzq(i, k);
          pt(i, k) = pow(ppres(i, k) * 1.e-5, rd_inv_cp) * pth(i, k);
          prho(i, k) = ppres(i, k) / (rd * pt(i, k));
          pinv_rho(i, k) = 1.0 / prho(i, k);
          prhofacr(i, k) = pow(rhosur * pinv_rho(i, k), 0.54);
        });
      team.team_barrier();

      bool log_qxpresent;
      const auto small = Globals<Real>::QSMALL;
      const int k_qxtop = find_top(team, Kokkos::subview(qr, i, Kokkos::ALL),
                                   small, kbot, ktop, kdir, log_qxpresent);

      if (log_qxpresent) {
        Real dt_left = dt;    // time remaining for sedi over full model (mp) time step
        Real prt_accum = 0.0; // precip rate for individual category

        int k_qxbot = find_bottom(team, Kokkos::subview(qr, i, Kokkos::ALL),
                                  small, kbot, k_qxtop, kdir, log_qxpresent);

        while (dt_left > 1.e-4) {
          Real Co_max = 0.0;
          int kmin, kmax;

          Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, m.num_vert / RealPack::n), [&] (int k) {
              pV_qr(i, k) = 0;
              pV_nr(i, k) = 0;
            });
          team.team_barrier();

          util::set_min_max(k_qxbot, k_qxtop, kmin, kmax);
          Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, kmax-kmin+1), [&] (int k_, Real& lmax) {
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
          team.team_barrier();

          // compute dt_sub
          const int tmpint1 = static_cast<int>(Co_max + 1.0);
          const Real dt_sub = util::min(dt_left, dt_left / tmpint1);

          // Move bottom cell down by 1 if not at ground already.
          const int k_temp = (k_qxbot == kbot) ? k_qxbot : k_qxbot - kdir;

          calc_first_order_upwind_step(m, team, i,
                                       k_temp, k_qxtop, kdir, dt_sub,
                                       sflux, sV_nr, snr);
          team.team_barrier();
          calc_first_order_upwind_step(m, team, i,
                                       k_temp, k_qxtop, kdir, dt_sub,
                                       sflux, sV_qr, sqr);
          team.team_barrier();

          // accumulated precip during time step
          if (k_qxbot == kbot) prt_accum += m.flux(i, kbot) * dt_sub;

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
