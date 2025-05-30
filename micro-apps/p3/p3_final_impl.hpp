#ifndef MICRO_SED_P3_FINAL_IMPL_HPP
#define MICRO_SED_P3_FINAL_IMPL_HPP

#include "p3_final.hpp"

#include "util.hpp"
#include "kokkos_util.hpp"
#include "wsm.hpp"
#include "micro_kokkos.hpp"
#include "p3_constants.hpp"
#include "scream_pack.hpp"
#include "p3_functions.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace p3 {
namespace micro_sed {

/*
 * Implementation of MicroSedFuncFinalKokkos, it implements things that were only
 * declared in p3_final.hpp. Clients should NOT #include this file, #include p3_final.hpp
 * instead.
 */

using scream::pack::IntPack;
using scream::pack::IntSmallPack;
using scream::pack::smallize;
using scream::pack::scalarize;
using scream::pack::BigPack;
using scream::pack::SmallPack;

template <typename ScalarT, typename DeviceT>
struct MicroSedFuncFinalKokkos<ScalarT,DeviceT>::Impl
{
  //
  // types
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;
  template <typename S> using SmallMask = scream::pack::Mask<SmallPack<S>::n>;

  using KT = KokkosTypes<Device>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  using view_1d_table = typename KT::template view_1d_table<Scalar, Globals<ScalarT>::MU_R_TABLE_DIM>;
  using view_2d_table = typename KT::template view_2d_table<Scalar, Globals<ScalarT>::VTABLE_DIM0,Globals<ScalarT>::VTABLE_DIM1>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  using ExeSpace    = typename KT::ExeSpace;
  using MemberType  = typename KT::MemberType;
  using TeamPolicy  = typename KT::TeamPolicy;

  using Fun = Functions<Scalar, Device>;

  //
  // members
  //

  static constexpr const char* NAME = "final";

private:
  int num_horz, num_vert, num_pack;

  view_2d_table vn_table, vm_table;
  view_1d_table mu_r_table;

  TeamPolicy policy;
  util::WorkspaceManager<Pack, Device> workspace_mgr;

public:
  Impl(int num_horz_, int num_vert_) :
    num_horz(num_horz_), num_vert(num_vert_),
    num_pack(scream::pack::npack<Pack>(num_vert_)),
    policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(num_horz, num_pack)),
    workspace_mgr(num_pack, 11, policy) // rain sed's high-water is 11 spaces for any team
  {
    Fun::init_kokkos_tables(vn_table, vm_table, mu_r_table);
  }

  int get_num_vert() const { return num_pack; }

  static void micro_sed_func (
    const Impl& msfk,
    const Int kts, const Int kte, const int its, const int ite, const Scalar dt,
    const view_2d<Pack>& qr, const view_2d<Pack>& nr,
    const view_2d<const Pack>& th, const view_2d<const Pack>& dzq, const view_2d<const Pack>& pres,
    const view_1d<Scalar>& prt_liq)
  {
    // constants
    const Scalar odt = 1.0 / dt;
    constexpr auto nsmall = Constants<Scalar>::NSMALL;
    constexpr auto rd = Constants<Scalar>::RD;
    constexpr auto rd_inv_cp = Constants<Scalar>::RD * Constants<Scalar>::INV_CP;
    constexpr auto rhosur = Constants<Scalar>::RHOSUR;
    constexpr auto qsmall = Constants<Scalar>::QSMALL;

    // direction of vertical leveling
    const Int kbot = (kts < kte) ? 0 : msfk.num_vert-1;
    const Int ktop = (kts < kte) ? msfk.num_vert-1 : 0;
    const Int kdir = (kts < kte) ? 1 : -1;

    const auto lqr = smallize(qr), lnr = smallize(nr);
    const auto sqr = scalarize(qr);

    // Rain sedimentation:  (adaptive substepping)
    // This parallel dispatch will give each vertical column to a thread team.
    // GPU architectures may have more threads than there are columns; in this
    // case, thread teams will contain multiple threads so that multiple threads
    // will be used to speed up inner (vertical over single-column) loops. The
    // TeamThreadRange dispatches below perform the intra-team parallelism.
    Kokkos::parallel_for(
      "main rain sed loop",
      msfk.policy,
      KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();

        auto workspace = msfk.workspace_mgr.get_workspace(team);

        // Get temporary workspaces needed for the rain-sed calculation
        Unmanaged<view_1d<Pack> > inv_dzq, rho, inv_rho, rhofacr, t, V_qr, V_nr, flux_qx, flux_nx, mu_r, lamr;
        workspace.template take_many_and_reset<11>(
          {"inv_dzq", "rho", "inv_rho", "rhofacr", "t", "V_qr", "V_nr", "flux_qx", "flux_nx", "mu_r", "lamr"},
          {&inv_dzq, &rho, &inv_rho, &rhofacr, &t, &V_qr, &V_nr, &flux_qx, &flux_nx, &mu_r, &lamr});

        // Get single-column subviews of all inputs, shouldn't need any i-indexing
        // after this.
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
        const Int k_qxtop = Fun::find_top(team, osqr, qsmall, kbot, ktop, kdir, log_qxpresent);

        if (log_qxpresent) {
          Scalar dt_left = dt;    // time remaining for sedi over full model (mp) time step
          Scalar prt_accum = 0.0; // precip rate for individual category

          Int k_qxbot = Fun::find_bottom(team, osqr, qsmall, kbot, k_qxtop, kdir, log_qxpresent);

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
                  typename Fun::Table3 table;
                  Spack tmp1, tmp2;
                  get_rain_dsd2(msfk.mu_r_table,
                                qr_gt_small, olqr(pk), olnr(pk), lmu_r(pk),
                                table.rdumii, table.dumii, llamr(pk),
                                tmp1, tmp2);
                  Fun::lookup(qr_gt_small, lmu_r(pk), llamr(pk), table);
                  // mass-weighted fall speed:
                  lV_qr(pk).set(qr_gt_small,
                                Fun::apply_table(qr_gt_small, msfk.vm_table, table) * lrhofacr(pk));
                  // number-weighted fall speed:
                  lV_nr(pk).set(qr_gt_small,
                                Fun::apply_table(qr_gt_small, msfk.vn_table, table) * lrhofacr(pk));
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

            Fun::template calc_first_order_upwind_step<2>(
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
              prt_liq(i) += prt_accum * Constants<Scalar>::INV_RHOW * odt;
            });
        }
      });

    // workspace_mgr.report(); // uncomment for detailed debug info
  }

  // Computes and returns rain size distribution parameters
  KOKKOS_INLINE_FUNCTION
  static void get_rain_dsd2 (
    const view_1d_table& mu_r_table,
    const SmallMask<Scalar>& qr_gt_small, const Spack& qr, Spack& nr, Spack& mu_r,
    Spack& rdumii, IntSmallPack& dumii, Spack& lamr,
    Spack& cdistr, Spack& logn0r)
  {
    constexpr auto nsmall = Constants<Scalar>::NSMALL;
    constexpr auto thrd = Constants<Scalar>::THRD;
    constexpr auto cons1 = Constants<Scalar>::CONS1;

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
          // Linearly interpolate mu_r.
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
};

template <typename Scalar, typename D>
MicroSedFuncFinalKokkos<Scalar, D>
::MicroSedFuncFinalKokkos (int num_horz, int num_vert) {
  impl = std::make_shared<Impl>(num_horz, num_vert);
}

template <typename Scalar, typename D>
void MicroSedFuncFinalKokkos<Scalar, D>
::micro_sed_func (
  const Int kts, const Int kte, const int its, const int ite, const Scalar dt,
  const view_2d<Pack>& qr, const view_2d<Pack>& nr,
  const view_2d<const Pack>& th, const view_2d<const Pack>& dzq, const view_2d<const Pack>& pres,
  const view_1d<Scalar>& prt_liq)
{
  Impl::micro_sed_func(*impl, kts, kte, its, ite, dt, qr, nr, th, dzq, pres,
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
