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

  KOKKOS_INLINE_FUNCTION
  static void lookup (
    const SmallMask<Scalar>& qr_gt_small, Table3& t, const SmallPack<Scalar>& mu_r, const SmallPack<Scalar>& lamr)
  {
    // find location in scaled mean size space
    const auto dum1 = (mu_r+1.) / lamr;
    const auto dum1_lt = qr_gt_small && (dum1 <= 195.e-6);
    if (dum1_lt.any()) {
      scream_masked_loop(dum1_lt, s) {
        const auto inv_dum3 = 0.1;
        auto rdumii = (dum1[s]*1.e6+5.)*inv_dum3;
        rdumii = util::max<Scalar>(rdumii,  1.);
        rdumii = util::min<Scalar>(rdumii, 20.);
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
      scream_masked_loop(dum1_gte, s) {
        const auto inv_dum3 = Globals<Scalar>::THRD*0.1;
        auto rdumii = (dum1[s]*1.e+6-195.)*inv_dum3 + 20.;
        rdumii = util::max<Scalar>(rdumii, 20.);
        rdumii = util::min<Scalar>(rdumii,300.);
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
  static Spack apply_table (
    const SmallMask<Scalar>& qr_gt_small, const view_2d_table& table, const Table3& t)
  {
    const auto rdumii_m_dumii = t.rdumii - Spack(t.dumii);
    const auto t_im1_jm1 = index(table, t.dumii-1, t.dumjj-1);
    const auto dum1 = (t_im1_jm1 + rdumii_m_dumii * t.inv_dum3 *
                       (index(table, t.dumii, t.dumjj-1) - t_im1_jm1));
    const auto t_im1_j = index(table, t.dumii-1, t.dumjj);
    const auto dum2 = (t_im1_j + rdumii_m_dumii * t.inv_dum3 *
                       (index(table, t.dumii, t.dumjj) - t_im1_j));
    return dum1 + (t.rdumjj - Spack(t.dumjj)) * (dum2 - dum1);
  }

  // Computes and returns rain size distribution parameters
  KOKKOS_INLINE_FUNCTION
  static void get_rain_dsd2_kokkos (
    const view_1d_table& mu_r_table,
    const SmallMask<Scalar>& qr_gt_small, const Spack& qr, Spack& nr, Spack& mu_r,
    Spack& rdumii, IntSmallPack& dumii, Spack& lamr,
    Spack& cdistr, Spack& logn0r)
  {
    constexpr auto nsmall = Globals<Scalar>::NSMALL;
    constexpr auto thrd = Globals<Scalar>::THRD;
    constexpr auto cons1 = Globals<Scalar>::CONS1;

    lamr = 0;
    cdistr = 0;
    logn0r = 0;

    // use lookup table to get mu
    // mu-lambda relationship is from Cao et al. (2008), eq. (7)

    // find spot in lookup table
    // (scaled N/q for lookup table parameter space)
    const auto nr_lim = max(nr, nsmall);
    Spack inv_dum(0);
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
        scream_masked_loop(m2, s) {
          // interpolate
          Scalar rdumiis = (inv_dum[s] - 250.e-6)*0.5e6;
          rdumiis = util::max<Scalar>(rdumiis, 1.0);
          rdumiis = util::min<Scalar>(rdumiis, 150.0);
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
      scream_masked_loop(either, s) {
        nr[s] = std::exp(3*std::log(lamr[s]) + std::log(qr[s]) +
                         std::log(std::tgamma(mu_r[s] + 1)) - std::log(std::tgamma(mu_r[s] + 4)))
          / cons1;
      }
    }
  }

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
    const view_1d_ptr_array<Spack, nfield>& r)
  {
    if (kdir == 1)
      calc_first_order_upwind_step< 1, nfield>(
        rho, inv_rho, inv_dzq, team, nk, k_bot, k_top, dt_sub, flux, V, r);
    else
      calc_first_order_upwind_step<-1, nfield>(
        rho, inv_rho, inv_dzq, team, nk, k_bot, k_top, dt_sub, flux, V, r);
  }

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
};

} // namespace micro_sed
} // namespace p3

#endif
