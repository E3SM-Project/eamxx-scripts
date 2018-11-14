#ifndef MICRO_SED_P3_FINAL_IMPL_HPP
#define MICRO_SED_P3_FINAL_IMPL_HPP

#include "p3_final.hpp"

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
  const SmallMask<Scalar>& qr_gt_small, Table3<Scalar>& t, const SmallPack<Scalar>& mu_r, const SmallPack<Scalar>& lamr)
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

template <typename Scalar, typename D>
struct MicroSedFuncFinalKokkos<Scalar,D>::Impl
{
  //
  // types
  //

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using KT = KokkosTypes<D>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  using view_1d_table = typename KT::template view_1d_table<Scalar, 150>;
  using view_2d_table = typename KT::template view_2d_table<Scalar, 300, 10>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  using ExeSpace    = typename KT::ExeSpace;
  using MemberType  = typename KT::MemberType;
  using TeamPolicy  = typename KT::TeamPolicy;

  //
  // members
  //

  static constexpr const char* NAME = "final";

private:
  int num_horz, num_vert, num_pack;

  view_2d_table vn_table, vm_table;
  view_1d_table mu_r_table;

  TeamPolicy policy;
  util::WorkspaceManager<Pack, D> workspace_mgr;

public:
  Impl(int num_horz_, int num_vert_) :
    num_horz(num_horz_), num_vert(num_vert_),
    num_pack(scream::pack::npack<Pack>(num_vert_)),
    vn_table("VN_TABLE"), vm_table("VM_TABLE"),
    mu_r_table("MU_R_TABLE"),
    policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(num_horz, num_pack)),
    workspace_mgr(num_pack, 11, policy) // rain sed's high-water is 11 workspace for any team
  {
    // initialize on host

    auto mirror_vn_table = Kokkos::create_mirror_view(vn_table);
    auto mirror_vm_table = Kokkos::create_mirror_view(vm_table);
    auto mirror_mu_table = Kokkos::create_mirror_view(mu_r_table);

    for (int i = 0; i < 300; ++i) {
      for (int k = 0; k < 10; ++k) {
        mirror_vn_table(i, k) = Globals<Scalar>::VN_TABLE[i][k];
        mirror_vm_table(i, k) = Globals<Scalar>::VM_TABLE[i][k];
      }
    }

    for (int i = 0; i < 150; ++i) {
      mirror_mu_table(i) = Globals<Scalar>::MU_R_TABLE[i];
    }

    // deep copy to device
    Kokkos::deep_copy(vn_table, mirror_vn_table);
    Kokkos::deep_copy(vm_table, mirror_vm_table);
    Kokkos::deep_copy(mu_r_table, mirror_mu_table);
  }

  int get_num_vert() const { return num_pack; }

  KOKKOS_INLINE_FUNCTION
  static Spack apply_table (
    const SmallMask<Scalar>& qr_gt_small, const view_2d_table& table, const Table3<Scalar>& t)
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

  static void micro_sed_func (
    const Impl& msfk,
    const Int kts, const Int kte, const int its, const int ite, const Scalar dt,
    const view_2d<Pack>& qr, const view_2d<Pack>& nr,
    const view_2d<Pack>& th, const view_2d<Pack>& dzq, const view_2d<Pack>& pres,
    const view_1d<Scalar>& prt_liq)
  {
    // constants
    const Scalar odt = 1.0 / dt;
    constexpr auto nsmall = Globals<Scalar>::NSMALL;
    constexpr auto rd = Globals<Scalar>::RD;
    constexpr auto rd_inv_cp = Globals<Scalar>::RD * Globals<Scalar>::INV_CP;
    constexpr auto rhosur = Globals<Scalar>::RHOSUR;
    constexpr auto qsmall = Globals<Scalar>::QSMALL;

    // direction of vertical leveling
    const Int kbot = (kts < kte) ? 0 : msfk.num_vert-1;
    const Int ktop = (kts < kte) ? msfk.num_vert-1 : 0;
    const Int kdir = (kts < kte) ? 1 : -1;

    const auto lqr = smallize(qr), lnr = smallize(nr);
    const auto sqr = scalarize(qr);

    // Rain sedimentation:  (adaptive substepping)
    Kokkos::parallel_for(
      "main rain sed loop",
      msfk.policy,
      KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();

        auto workspace = msfk.workspace_mgr.get_workspace(team);

        Unmanaged<view_1d<Pack> > inv_dzq, rho, inv_rho, rhofacr, t, V_qr, V_nr, flux_qx, flux_nx, mu_r, lamr;
        workspace.template take_many_and_reset<11>(
          {"inv_dzq", "rho", "inv_rho", "rhofacr", "t", "V_qr", "V_nr", "flux_qx", "flux_nx", "mu_r", "lamr"},
          {&inv_dzq, &rho, &inv_rho, &rhofacr, &t, &V_qr, &V_nr, &flux_qx, &flux_nx, &mu_r, &lamr});

        const auto olqr = util::subview(lqr, i), olnr = util::subview(lnr, i);
        const auto osqr = util::subview(sqr, i);
        const auto odzq = util::subview(dzq, i), oth = util::subview(th, i), opres = util::subview(pres, i);
        const auto
          lflux_qx = smallize(flux_qx), lflux_nx = smallize(flux_nx),
          lmu_r = smallize(mu_r), llamr = smallize(lamr),
          lrho = smallize(rho), linv_rho = smallize(inv_rho), lrhofacr = smallize(rhofacr),
          linv_dzq = smallize(inv_dzq);

        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, msfk.num_pack), [&] (Int k) {
            inv_dzq(k) = 1.0 / odzq(k);
            t(k) = pow(opres(k) * 1.e-5, rd_inv_cp) * oth(k);
            rho(k) = opres(k) / (rd * t(k));
            inv_rho(k) = 1.0 / rho(k);
            rhofacr(k) = pow(rhosur * inv_rho(k), 0.54);
          });
        team.team_barrier();

        bool log_qxpresent;
        const Int k_qxtop = find_top(team, osqr, qsmall, kbot, ktop, kdir, log_qxpresent);

        if (log_qxpresent) {
          Scalar dt_left = dt;    // time remaining for sedi over full model (mp) time step
          Scalar prt_accum = 0.0; // precip rate for individual category

          Int k_qxbot = find_bottom(team, osqr, qsmall, kbot, k_qxtop, kdir, log_qxpresent);

          const auto lV_qr = smallize(V_qr), lV_nr = smallize(V_nr);
          const auto sflux_qx = scalarize(lflux_qx);

          while (dt_left > 1.e-4) {
            Scalar Co_max = 0.0;
            Int kmin, kmax;

            Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, msfk.num_pack), [&] (Int k) {
                V_qr(k) = 0;
                V_nr(k) = 0;
              });
            team.team_barrier();

            util::set_min_max(k_qxbot, k_qxtop, kmin, kmax, Spack::n);

            Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(team, kmax-kmin+1), [&] (int pk_, Scalar& lmax) {
                const int pk = kmin + pk_;
                auto qr_gt_small = (olqr(pk) > qsmall);
                if (qr_gt_small.any()) {
                  // Compute Vq, Vn:
                  olnr(pk).set(qr_gt_small, max(olnr(pk), nsmall));
                  Table3<Scalar> table;
                  Spack tmp1, tmp2;
                  get_rain_dsd2_kokkos(msfk.mu_r_table,
                                       qr_gt_small, olqr(pk), olnr(pk), lmu_r(pk),
                                       table.rdumii, table.dumii, llamr(pk),
                                       tmp1, tmp2);
                  find_lookupTable_indices_3_kokkos(qr_gt_small, table, lmu_r(pk), llamr(pk));
                  // mass-weighted fall speed:
                  lV_qr(pk).set(qr_gt_small,
                                apply_table(qr_gt_small, msfk.vm_table, table) * lrhofacr(pk));
                  // number-weighted fall speed:
                  lV_nr(pk).set(qr_gt_small,
                                apply_table(qr_gt_small, msfk.vn_table, table) * lrhofacr(pk));
                  const auto Co_max_local = max(qr_gt_small, -1,
                                                lV_qr(pk) * dt_left * linv_dzq(pk));
                  if (Co_max_local > lmax)
                    lmax = Co_max_local;
                }
              }, Kokkos::Max<Scalar>(Co_max));
            team.team_barrier();

            if (Co_max < 0) {
              // qr is everywhere too small. Exit dt_left loop.
              break;
            }

            // compute dt_sub
            const Int Co_max_p1 = static_cast<Int>(Co_max + 1.0);
            const Scalar dt_sub = util::min(dt_left, dt_left / Co_max_p1);

            // Move bottom cell down by 1 if not at ground already.
            const Int k_temp = (k_qxbot == kbot) ? k_qxbot : k_qxbot - kdir;

            calc_first_order_upwind_step<2>(
              lrho, linv_rho, linv_dzq, team,
              msfk.num_vert, k_temp, k_qxtop, kdir, dt_sub,
              {&lflux_qx, &lflux_nx}, {&lV_qr, &lV_nr}, {&olqr, &olnr});
            team.team_barrier();

            // accumulated precip during time step
            if (k_qxbot == kbot) prt_accum += sflux_qx(kbot) * dt_sub;

            dt_left -= dt_sub;  // update time remaining for sedimentation
            if (k_qxbot != kbot) k_qxbot -= kdir;
          }

          Kokkos::single(
            Kokkos::PerTeam(team), [&] () {
              prt_liq(i) += prt_accum * Globals<Scalar>::INV_RHOW * odt;
            });
        }
      });

    // workspace_mgr.report(); // uncomment for detailed debug info
  }

};

template <typename Scalar, typename D>
MicroSedFuncFinalKokkos<Scalar, D>
::MicroSedFuncFinalKokkos (int num_horz, int num_vert) {
  impl = std::make_shared<Impl>(num_horz, num_vert);
}

template <typename Scalar, typename D>
void MicroSedFuncFinalKokkos<Scalar, D>
::micro_sed_func (
  MicroSedFuncFinalKokkos& msfk,
  const Int kts, const Int kte, const int its, const int ite, const Scalar dt,
  const view_2d<Pack>& qr, const view_2d<Pack>& nr,
  const view_2d<Pack>& th, const view_2d<Pack>& dzq, const view_2d<Pack>& pres,
  const view_1d<Scalar>& prt_liq)
{
  Impl::micro_sed_func(*msfk.impl, kts, kte, its, ite, dt, qr, nr, th, dzq, pres,
                       prt_liq);
}

template <typename Scalar, typename D>
int MicroSedFuncFinalKokkos<Scalar, D>
::get_num_vert () const {
  return impl->get_num_vert();
}

} // namespace p3
} // namespace micro_sed

#endif
