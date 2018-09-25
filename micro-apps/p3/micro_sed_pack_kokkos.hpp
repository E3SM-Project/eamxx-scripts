#ifndef MICRO_SED_PACK_KOKKOS_HPP
#define MICRO_SED_PACK_KOKKOS_HPP

#include "util.hpp"
#include "initial_conditions.hpp"
#include "micro_kokkos.hpp"
#include "micro_sed_vanilla.hpp"
#include "micro_sed_vanilla_kokkos.hpp"
#include "scream_pack.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace p3 {
namespace micro_sed {

using micro_sed::Globals;

template <typename T>
using BigPack = scream::pack::Pack<T, SCREAM_PACKN>;
template <typename T>
using SmallPack = scream::pack::Pack<T, SCREAM_PACKN / SCREAM_SMALL_PACK_FACTOR>;
using Mask = scream::pack::Mask<BigPack<Real>::n>;
using SmallMask = scream::pack::Mask<SmallPack<Real>::n>;

using RealPack = BigPack<Real>;
using IntPack = BigPack<Int>;
using RealSmallPack = SmallPack<Real>;
using IntSmallPack = SmallPack<Int>;

template <typename T, typename ...Parms> KOKKOS_FORCEINLINE_FUNCTION
Kokkos::View<T**, Parms..., Kokkos::MemoryTraits<Kokkos::Unmanaged> >
scalarize (const Kokkos::View<BigPack<T>**, Parms...>& vp) {
  return Kokkos::View<T**, Parms..., Kokkos::MemoryTraits<Kokkos::Unmanaged> >(
    reinterpret_cast<T*>(vp.data()), vp.extent_int(0), RealPack::n * vp.extent_int(1));
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
struct MicroSedFuncPackKokkos : public MicroSedFuncVanillaKokkos<Real, kokkos_2d_t<RealPack> > {
  int num_pack;

  static constexpr const char* NAME = "kokkos_pack";

public:
  MicroSedFuncPackKokkos(int num_horz_, int num_vert_) :
    MicroSedFuncVanillaKokkos<Real, kokkos_2d_t<RealPack> >(num_horz_, scream::pack::npack<RealPack>(num_vert_)),
    num_pack(scream::pack::npack<RealPack>(num_vert_))
  {
    this->num_vert = num_vert_;
  }

  static std::string custom_msg()
  {
    std::ostringstream out;
    out << " " << SCREAM_PACKN << " " << SCREAM_SMALL_PACK_FACTOR;
    return out.str();
  }

  virtual int get_num_vert() const { return num_pack; }
};

struct Table3 {
  IntSmallPack dumii, dumjj;
  RealSmallPack rdumii, rdumjj, inv_dum3;
};

KOKKOS_INLINE_FUNCTION
void find_lookupTable_indices_3_kokkos (
  const SmallMask& qr_gt_small, Table3& t, const RealSmallPack& mu_r, const RealSmallPack& lamr)
{
  // find location in scaled mean size space
  const auto dum1 = (mu_r+1.) / lamr;
  const auto dum1_lt = qr_gt_small && (dum1 <= 195.e-6);
  if (dum1_lt.any()) {
    scream_masked_loop(dum1_lt) {
      const auto inv_dum3 = 0.1;
      auto rdumii = (dum1[s]*1.e6+5.)*inv_dum3;
      rdumii = util::max<Real>(rdumii,  1.);
      rdumii = util::min<Real>(rdumii, 20.);
      Int dumii = rdumii;
      dumii = util::max(dumii,  1);
      dumii = util::min(dumii, 20);
      t.inv_dum3[s] = inv_dum3;
      t.rdumii[s] = rdumii;
      t.dumii[s] = dumii;
    }
  }
  const auto dum1_gte = qr_gt_small && ! dum1_lt;
  if (dum1_gte.any()) {
    scream_masked_loop(dum1_gte) {
      const auto inv_dum3 = Globals<Real>::THRD*0.1;
      auto rdumii = (dum1[s]*1.e+6-195.)*inv_dum3 + 20.;
      rdumii = util::max<Real>(rdumii, 20.);
      rdumii = util::min<Real>(rdumii,300.);
      Int dumii = rdumii;
      dumii = util::max(dumii, 20);
      dumii = util::min(dumii,299);
      t.inv_dum3[s] = inv_dum3;
      t.rdumii[s] = rdumii;
      t.dumii[s] = dumii;
    }
  }

  // find location in mu_r space
  {
    auto rdumjj = mu_r+1.;
    rdumjj = max(rdumjj,1.);
    rdumjj = min(rdumjj,10.);
    IntSmallPack dumjj(rdumjj);
    dumjj  = max(dumjj,1);
    dumjj  = min(dumjj,9);
    t.rdumjj.set(qr_gt_small, rdumjj);
    t.dumjj.set(qr_gt_small, dumjj);
  }
}

KOKKOS_INLINE_FUNCTION
RealSmallPack apply_table (
  const SmallMask& qr_gt_small, const kokkos_2d_table_t<Real>& table, const Table3& t)
{
  const auto rdumii_m_dumii = t.rdumii-RealSmallPack(t.dumii);
  const auto t_im1_jm1 = index(table, t.dumii-1, t.dumjj-1);
  const auto dum1 = (t_im1_jm1 + rdumii_m_dumii * t.inv_dum3 *
                     (index(table, t.dumii, t.dumjj-1) - t_im1_jm1));
  const auto t_im1_j = index(table, t.dumii-1, t.dumjj);
  const auto dum2 = (t_im1_j + rdumii_m_dumii * t.inv_dum3 *
                     (index(table, t.dumii, t.dumjj) - t_im1_j));
  return dum1 + (t.rdumjj - RealSmallPack(t.dumjj)) * (dum2 - dum1);
}

// Computes and returns rain size distribution parameters
KOKKOS_INLINE_FUNCTION
void get_rain_dsd2_kokkos (
  const SmallMask& qr_gt_small, const RealSmallPack& qr, RealSmallPack& nr, RealSmallPack& mu_r,
  RealSmallPack& rdumii, IntSmallPack& dumii, RealSmallPack& lamr,
  const kokkos_1d_table_t<Real>& mu_r_table, RealSmallPack& cdistr, RealSmallPack& logn0r)
{
  constexpr auto nsmall = Globals<Real>::NSMALL;
  constexpr auto thrd = Globals<Real>::THRD;
  constexpr auto cons1 = Globals<Real>::CONS1;

  lamr = 0;
  cdistr = 0;
  logn0r = 0;

  // use lookup table to get mu
  // mu-lambda relationship is from Cao et al. (2008), eq. (7)

  // find spot in lookup table
  // (scaled N/q for lookup table parameter space)
  const auto nr_lim = max(nr, nsmall);
  RealSmallPack inv_dum(0);
  inv_dum.set(qr_gt_small,
              pow(qr / (cons1 * nr_lim * 6.0), thrd));

  mu_r = 0;
  {
    const auto m1 = qr_gt_small && (inv_dum < 282.e-6);
    mu_r.set(m1, 8.282);
  }
  {
    const auto m2 = qr_gt_small && (inv_dum >= 282.e-6) && (inv_dum < 502.e-6);
    if (m2.any()) {
      scream_masked_loop(m2) {
        // interpolate
        Real rdumiis = (inv_dum[s] - 250.e-6)*0.5e6;
        rdumiis = util::max<Real>(rdumiis, 1.0);
        rdumiis = util::min<Real>(rdumiis, 150.0);
        rdumii[s] = rdumiis;
        Int dumiis = rdumiis;
        dumiis = util::min(dumiis, 149);
        dumii[s] = dumiis;
        const auto mu_r_im1 = mu_r_table(dumiis-1);
        mu_r[s] = mu_r_im1 + (mu_r_table(dumiis) - mu_r_im1) * (rdumiis - dumiis);
      }
    }
  }

  // recalculate slope based on mu_r
  lamr.set(qr_gt_small,
           pow(cons1 * nr_lim * (mu_r + 3) *
               (mu_r + 2) * (mu_r + 1)/qr,
               thrd));

  // check for slope
  const auto lammax = (mu_r+1.)*1.e+5;
  // set to small value since breakup is explicitly included (mean size 0.8 mm)
  const auto lammin = (mu_r+1.)*1250.0;
  // apply lambda limiters for rain
  const auto lt = qr_gt_small && (lamr < lammin);
  const auto gt = qr_gt_small && (lamr > lammax);
  const auto either = lt || gt;
  nr.set(qr_gt_small, nr_lim);
  if (either.any()) {
    lamr.set(lt, lammin);
    lamr.set(gt, lammax);
    scream_masked_loop(either) {
      nr[s] = std::exp(3*std::log(lamr[s]) + std::log(qr[s]) +
                       std::log(std::tgamma(mu_r[s] + 1)) - std::log(std::tgamma(mu_r[s] + 4)))
        / cons1;
    }
  }
}

//TODO Unit test.
// Calculate the step in the region [k_bot, k_top].
template <Int kdir, int nfield, typename MSPK>
KOKKOS_INLINE_FUNCTION
void calc_first_order_upwind_step (
  const MSPK& m, const member_type& team, const int& i,
  const Int& nk, const Int& k_bot, const Int& k_top, const Real& dt_sub,
  const Kokkos::Array<const kokkos_2d_t<RealSmallPack>*,nfield>& flux,
  const Kokkos::Array<const kokkos_2d_t<RealSmallPack>*,nfield>& V,
  const Kokkos::Array<const kokkos_2d_t<RealSmallPack>*,nfield>& r)
{
  const kokkos_2d_t<RealSmallPack>
    rho = smallize(m.rho),
    inv_rho = smallize(m.inv_rho),
    inv_dzq = smallize(m.inv_dzq);

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
        (*flux[f])(i, k) = (*V[f])(i, k) * (*r[f])(i, k) * rho(i, k);
    });
  team.team_barrier();

  Kokkos::single(
    Kokkos::PerTeam(team), [&] () {
      const Int k = k_top_pack;
      if (nk % RealSmallPack::n != 0) {
        const auto mask =
          scream::pack::range<IntSmallPack>(k_top_pack*RealSmallPack::n) >= nk;
        for (int f = 0; f < nfield; ++f)
          (*flux[f])(i, k_top_pack).set(mask, 0);
      }
      for (int f = 0; f < nfield; ++f) {
        // compute flux divergence
        const auto flux_pkdir = (kdir == -1) ?
          shift_right(0, (*flux[f])(i, k)) :
          shift_left (0, (*flux[f])(i, k));
        const auto fluxdiv = (flux_pkdir - (*flux[f])(i, k)) * inv_dzq(i, k);
        // update prognostic variables
        (*r[f])(i, k) += fluxdiv * dt_sub * inv_rho(i, k);
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
          shift_right((*flux[f])(i, k+kdir), (*flux[f])(i, k)) :
          shift_left ((*flux[f])(i, k+kdir), (*flux[f])(i, k));
        const auto fluxdiv = (flux_pkdir - (*flux[f])(i, k)) * inv_dzq(i, k);
        // update prognostic variables
        (*r[f])(i, k) += fluxdiv * dt_sub * inv_rho(i, k);
      }
    });
}

template <int nfield, typename MSPK> KOKKOS_INLINE_FUNCTION
void calc_first_order_upwind_step (
  const MSPK& m, const member_type& team, const int& i,
  const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, const Real& dt_sub,
  const Kokkos::Array<const kokkos_2d_t<RealSmallPack>*,nfield>& flux,
  const Kokkos::Array<const kokkos_2d_t<RealSmallPack>*,nfield>& V,
  const Kokkos::Array<const kokkos_2d_t<RealSmallPack>*,nfield>& r)
{
  if (kdir == 1)
    calc_first_order_upwind_step< 1, nfield>(
      m, team, i, nk, k_bot, k_top, dt_sub, flux, V, r);
  else
    calc_first_order_upwind_step<-1, nfield>(
      m, team, i, nk, k_bot, k_top, dt_sub, flux, V, r);
}

//TODO Unit test.
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

template <typename MSPK>
void micro_sed_func (
  MSPK& m,
  const Int kts, const Int kte, const int its, const int ite, const Real dt,
  const kokkos_2d_t<RealPack>& qr, const kokkos_2d_t<RealPack>& nr,
  const kokkos_2d_t<RealPack>& th, const kokkos_2d_t<RealPack>& dzq, const kokkos_2d_t<RealPack>& pres,
  const kokkos_1d_t<Real>& prt_liq)
{
  const kokkos_2d_t<RealSmallPack>
    linv_dzq = smallize(m.inv_dzq),
    lrhofacr = smallize(m.rhofacr),
    lV_qr = smallize(m.V_qr),
    lV_nr = smallize(m.V_nr),
    lflux_qx = smallize(m.flux_qx),
    lflux_nx = smallize(m.flux_nx),
    lqr = smallize(qr),
    lnr = smallize(nr),
    lmu_r = smallize(m.mu_r),
    llamr = smallize(m.lamr);
  const kokkos_2d_t<Real>
    sqr = scalarize(qr),
    sflux_qx = scalarize(m.flux_qx);

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
    util::ExeSpaceUtils<>::get_default_team_policy(m.num_horz, m.num_pack),
    KOKKOS_LAMBDA(const member_type& team) {
      const int i = team.league_rank();

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, m.num_pack), [&] (Int k) {
          m.inv_dzq(i, k) = 1.0 / dzq(i, k);
          m.t(i, k) = pow(pres(i, k) * 1.e-5, rd_inv_cp) * th(i, k);
          m.rho(i, k) = pres(i, k) / (rd * m.t(i, k));
          m.inv_rho(i, k) = 1.0 / m.rho(i, k);
          m.rhofacr(i, k) = pow(rhosur * m.inv_rho(i, k), 0.54);
        });
      team.team_barrier();

      bool log_qxpresent;
      const Int k_qxtop = find_top(team, Kokkos::subview(sqr, i, Kokkos::ALL),
                                   qsmall, kbot, ktop, kdir, log_qxpresent);

      if (log_qxpresent) {
        Real dt_left = dt;    // time remaining for sedi over full model (mp) time step
        Real prt_accum = 0.0; // precip rate for individual category

        Int k_qxbot = find_bottom(team, Kokkos::subview(sqr, i, Kokkos::ALL),
                                  qsmall, kbot, k_qxtop, kdir, log_qxpresent);

        while (dt_left > 1.e-4) {
          Real Co_max = 0.0;
          Int kmin, kmax;

          Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, m.num_pack), [&] (Int k) {
              m.V_qr(i, k) = 0;
              m.V_nr(i, k) = 0;
            });
          team.team_barrier();

          util::set_min_max(k_qxbot, k_qxtop, kmin, kmax, RealSmallPack::n);
          Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, kmax-kmin+1), [&] (int pk_, Real& lmax) {
              const int pk = kmin + pk_;
              auto qr_gt_small = (lqr(i, pk) > qsmall);
              if (qr_gt_small.any()) {
                // Compute Vq, Vn:
                lnr(i, pk).set(qr_gt_small, max(lnr(i, pk), nsmall));
                Table3 t;
                RealSmallPack tmp1, tmp2;
                get_rain_dsd2_kokkos(qr_gt_small, lqr(i, pk), lnr(i, pk), lmu_r(i, pk),
                                     t.rdumii, t.dumii, llamr(i, pk),
                                     m.mu_r_table, tmp1, tmp2);
                find_lookupTable_indices_3_kokkos(qr_gt_small, t, lmu_r(i, pk), llamr(i, pk));
                // mass-weighted fall speed:
                lV_qr(i, pk).set(qr_gt_small,
                                 apply_table(qr_gt_small, m.vm_table, t) * lrhofacr(i, pk));
                // number-weighted fall speed:
                lV_nr(i, pk).set(qr_gt_small,
                                 apply_table(qr_gt_small, m.vn_table, t) * lrhofacr(i, pk));
                const auto Co_max_local = max(qr_gt_small, -1,
                                              lV_qr(i, pk) * dt_left * linv_dzq(i, pk));
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
            m, team, i,
            m.num_vert, k_temp, k_qxtop, kdir, dt_sub,
            {&lflux_qx, &lflux_nx}, {&lV_qr, &lV_nr}, {&lqr, &lnr});
          team.team_barrier();

          // accumulated precip during time step
          if (k_qxbot == kbot) prt_accum += sflux_qx(i, kbot) * dt_sub;

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

template <typename Real>
void populate_kokkos_from_vec (
  const int num_horz, const int num_vert, vector_2d_t<Real> const& vec, kokkos_2d_t<RealPack>& device)
{
  const auto mirror = Kokkos::create_mirror_view(device);
  const auto smirror = scalarize(mirror);

  for (int i = 0; i < num_horz; ++i) {
    for (Int k = 0; k < num_vert; ++k) {
      smirror(i, k) = vec[i][k];
    }
  }

  Kokkos::deep_copy(device, mirror);
}

template <typename Real>
void dump_to_file_k (
  const char* basename,
  const kokkos_2d_t<RealPack>& qr, const kokkos_2d_t<RealPack>& nr, const kokkos_2d_t<RealPack>& th,
  const kokkos_2d_t<RealPack>& dzq, const kokkos_2d_t<RealPack>& pres, const kokkos_1d_t<Real>& prt_liq,
  const int ni, const int nk, const Real dt, const int ts)
{
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

  const auto
    sqr = scalarize(qr_m),
    snr = scalarize(nr_m),
    sth = scalarize(th_m),
    sdzq = scalarize(dzq_m),
    spres = scalarize(pres_m);
  const int ldk = sqr.extent_int(1);

  util::dump_to_file(
    basename, sqr.data(), snr.data(), sth.data(), sdzq.data(), spres.data(),
    prt_liq_m.data(), ni, nk, dt, ts, ldk);
}

} // namespace p3
} // namespace micro_sed

#endif
