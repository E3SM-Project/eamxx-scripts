#ifndef MICRO_SED_VANILLA_KOKKOS_HPP
#define MICRO_SED_VANILLA_KOKKOS_HPP

#include "util.hpp"
#include "initial_conditions.hpp"
#include "micro_kokkos.hpp"
#include "p3_common.hpp"

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>

namespace p3 {
namespace micro_sed {

/**
 * Finds indices in rain lookup table (3)
 */
template <typename Scalar> KOKKOS_FUNCTION
void find_lookupTable_indices_3_kokkos(int& dumii, int& dumjj, Scalar& rdumii, Scalar& rdumjj, Scalar& inv_dum3,
                                       const Scalar mu_r, const Scalar lamr)
{
  // find location in scaled mean size space
  Scalar dum1 = (mu_r+1.) / lamr;
  if (dum1 <= 195.e-6) {
    inv_dum3  = 0.1;
    rdumii = (dum1*1.e6+5.)*inv_dum3;
    rdumii = util::max<Scalar>(rdumii, 1.);
    rdumii = util::min<Scalar>(rdumii,20.);
    dumii  = static_cast<int>(rdumii);
    dumii  = util::max(dumii, 1);
    dumii  = util::min(dumii,20);
  }
  else {
    inv_dum3  = Globals<Scalar>::THRD*0.1;           // i.e. 1/30
    rdumii = (dum1*1.e+6-195.)*inv_dum3 + 20.;
    rdumii = util::max<Scalar>(rdumii, 20.);
    rdumii = util::min<Scalar>(rdumii,300.);
    dumii  = static_cast<int>(rdumii);
    dumii  = util::max(dumii, 20);
    dumii  = util::min(dumii,299);
  }

  // find location in mu_r space
  rdumjj = mu_r+1.;
  rdumjj = util::max<Scalar>(rdumjj,1.);
  rdumjj = util::min<Scalar>(rdumjj,10.);
  dumjj  = static_cast<int>(rdumjj);
  dumjj  = util::max(dumjj,1);
  dumjj  = util::min(dumjj,9);
}

template <typename Scalar, typename D=DefaultDevice>
struct MicroSedFuncVanillaKokkos
{
  //
  // types
  //

  using Pack = Scalar;

  template <typename S>
  using view_1d = typename KokkosTypes<D>::template view_1d<S>;
  template <typename S>
  using view_2d = typename KokkosTypes<D>::template view_2d<S>;

  using ExeSpace    = typename KokkosTypes<D>::ExeSpace;
  using MemberType  = typename KokkosTypes<D>::MemberType;

  using view_1d_table = typename KokkosTypes<D>::template view_1d_table<Scalar, 150>;
  using view_2d_table = typename KokkosTypes<D>::template view_2d_table<Scalar, 300, 10>;

  //
  // members
  //

private:
  int num_horz, num_vert;

  // re-usable scratch views
  view_2d<Scalar> V_qr, V_nr, flux_qx, flux_nx, mu_r, lamr, rhofacr, inv_dzq, rho, inv_rho, t;

  view_2d_table vn_table, vm_table;
  view_1d_table mu_r_table;

public:
  static constexpr const char* NAME = "vanilla";

  MicroSedFuncVanillaKokkos(int num_horz_, int num_vert_) :
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
    vn_table("VN_TABLE"), vm_table("VM_TABLE"),
    mu_r_table("MU_R_TABLE")
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

  static std::string custom_msg()
  {
    return "";
  }

  int get_num_vert() const { return num_vert; }

  static void reset(const MicroSedFuncVanillaKokkos<Scalar, D>& msvk)
  {
    Kokkos::parallel_for(
      "2d reset",
      util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(msvk.num_horz, msvk.num_vert),
      KOKKOS_LAMBDA(MemberType team_member) {
        const int i = team_member.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int k) {
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
        });
      });
  }

  /**
   * Computes and returns rain size distribution parameters
   */
  KOKKOS_FUNCTION
  static void get_rain_dsd2_kokkos(
    const view_1d_table& mu_r_table,
    const Scalar qr, Scalar& nr, Scalar& mu_r, Scalar& rdumii, int& dumii, Scalar& lamr,
    Scalar& cdistr, Scalar& logn0r)
  {
    constexpr Scalar nsmall = Globals<Scalar>::NSMALL;
    if (qr >= Globals<Scalar>::QSMALL) {
      // use lookup table to get mu
      // mu-lambda relationship is from Cao et al. (2008), eq. (7)

      // find spot in lookup table
      // (scaled N/q for lookup table parameter space_
      nr = util::max(nr, nsmall);
      Scalar inv_dum = std::pow(qr / (Globals<Scalar>::CONS1 * nr * 6.0), Globals<Scalar>::THRD);

      if (inv_dum < 282.e-6) {
        mu_r = 8.282;
      }
      else if (inv_dum >= 282.e-6 && inv_dum < 502.e-6) {
        // interpolate
        rdumii = (inv_dum-250.e-6)*1.e+6*0.5;
        rdumii = util::max<Scalar>(rdumii,1.0);
        rdumii = util::min<Scalar>(rdumii,150.0);
        dumii  = static_cast<int>(rdumii);
        dumii  = util::min(149,dumii);
        mu_r   = mu_r_table(dumii-1) + (mu_r_table(dumii) - mu_r_table(dumii-1)) * (rdumii-dumii);
      }
      else if (inv_dum >= 502.e-6) {
        mu_r = 0.0;
      }

      lamr   = std::pow((Globals<Scalar>::CONS1 *nr *(mu_r+3.0) * (mu_r+2) * (mu_r+1.)/(qr)), Globals<Scalar>::THRD); // recalculate slope based on mu_r
      Scalar lammax = (mu_r+1.)*1.e+5;  // check for slope
      Scalar lammin = (mu_r+1.)*1250.0; // set to small value since breakup is explicitly included (mean size 0.8 mm)

      // apply lambda limiters for rain
      if (lamr < lammin) {
        lamr = lammin;
        nr   = std::exp(3.*std::log(lamr) + std::log(qr) + std::log(std::tgamma(mu_r+1.)) - std::log(std::tgamma(mu_r+4.)))/(Globals<Scalar>::CONS1);
      }
      else if (lamr > lammax) {
        lamr = lammax;
        nr   = std::exp(3.*std::log(lamr) + std::log(qr) + std::log(std::tgamma(mu_r+1.)) - log(std::tgamma(mu_r+4.)))/(Globals<Scalar>::CONS1);
      }

      cdistr  = 0; //nr/std::tgamma(mu_r+1.);
      logn0r  = 0; //std::log10(nr) + (mu_r+1.)*std::log10(lamr) - std::log10(std::tgamma(mu_r+1)); // note: logn0r is calculated as log10(n0r);
    }
    else {
      lamr   = 0.0;
      cdistr = 0.0;
      logn0r = 0.0;
    }
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
  void micro_sed_func(
    const int kts, const int kte, const int its, const int ite, const Scalar dt,
    view_2d<Scalar>& qr, view_2d<Scalar>& nr,
    view_2d<const Scalar> const& th, view_2d<const Scalar> const& dzq, view_2d<const Scalar> const& pres,
    view_1d<Scalar>& prt_liq)
  {
    micro_sed_func_impl(*this, kts, kte, its, ite, dt, qr, nr, th, dzq, pres, prt_liq);
  }

  static void micro_sed_func_impl(
    const MicroSedFuncVanillaKokkos<Scalar, D>& msvk,
    const int kts, const int kte, const int its, const int ite, const Scalar dt,
    view_2d<Scalar>& qr, view_2d<Scalar>& nr,
    view_2d<const Scalar> const& th, view_2d<const Scalar> const& dzq, view_2d<const Scalar> const& pres,
    view_1d<Scalar>& prt_liq)
  {
    // constants
    const Scalar odt = 1.0 / dt;
    constexpr Scalar nsmall = Globals<Scalar>::NSMALL;

    // direction of vertical leveling
    ConstExceptGnu int kbot = (kts < kte) ? 0: msvk.num_vert-1;
    ConstExceptGnu int kdir = (kts < kte) ? 1  : -1;

    // Rain sedimentation:  (adaptivive substepping)
    Kokkos::parallel_for(
      "main rain sed loop",
      util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(msvk.num_horz, msvk.num_vert),
      KOKKOS_LAMBDA(MemberType team_member) {
        const int i = team_member.league_rank();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int k) {
          // inverse of thickness of layers
          msvk.inv_dzq(i, k) = 1 / dzq(i, k);
          msvk.t(i, k) = std::pow(pres(i, k) * 1.e-5, Globals<Scalar>::RD * Globals<Scalar>::INV_CP) * th(i, k);
          msvk.rho(i, k) = pres(i, k) / (Globals<Scalar>::RD * msvk.t(i, k));
          msvk.inv_rho(i, k) = 1.0 / msvk.rho(i, k);
          msvk.rhofacr(i, k) = std::pow(Globals<Scalar>::RHOSUR * msvk.inv_rho(i, k), 0.54);
        });
        team_member.team_barrier();

        // Note, we are skipping supersaturation checks

        bool log_qxpresent = false;
        int k_qxtop = -1; // avoid warning, but don't use a meanigful value

        // find top, determine qxpresent
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int k, int& lmax) {
          if (qr(i, k) >= Globals<Scalar>::QSMALL && k*kdir > lmax) {
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

          Scalar dt_left = dt;    // time remaining for sedi over full model (mp) time step
          Scalar prt_accum = 0.0; // precip rate for individual category
          int k_qxbot = 0;

          // find bottom
          Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int k, int& lmin) {
            if (qr(i, k) >= Globals<Scalar>::QSMALL && k*kdir < lmin) {
              lmin = k*kdir;
            }
          }, Kokkos::Min<int>(k_qxbot));
          // As log_qxpresent is true, we don't have to worry about this
          // reduction as we did for the one for k_qxtop.
          k_qxbot *= kdir;

          while (dt_left > 1.e-4) {
            Scalar Co_max = 0.0;
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, msvk.num_vert), [=] (int kk) {
              msvk.V_qr(i, kk) = 0.0;
              msvk.V_nr(i, kk) = 0.0;
            });
            team_member.team_barrier();

            int kmin, kmax;
            util::set_min_max(k_qxtop, k_qxbot, kmin, kmax);
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, kmax-kmin+1), [=] (int k_, Scalar& lmax) {
              const int k = kmin + k_;
              if (qr(i, k) > Globals<Scalar>::QSMALL) {
                // Compute Vq, Vn:
                nr(i, k) = util::max(nr(i, k), nsmall);
                Scalar rdumii=0.0, tmp1=0.0, tmp2=0.0, rdumjj=0.0, inv_dum3=0.0;
                int dumii=0, dumjj=0;
                get_rain_dsd2_kokkos(msvk.mu_r_table, qr(i, k), nr(i, k), msvk.mu_r(i, k), rdumii, dumii, msvk.lamr(i, k), tmp1, tmp2);
                find_lookupTable_indices_3_kokkos(dumii, dumjj, rdumii, rdumjj, inv_dum3, msvk.mu_r(i, k), msvk.lamr(i, k));

                // mass-weighted fall speed:
                Scalar dum1 = msvk.vm_table(dumii-1, dumjj-1) + (rdumii-dumii) * inv_dum3 *
                  (msvk.vm_table(dumii, dumjj-1) - msvk.vm_table(dumii-1, dumjj-1));
                Scalar dum2 = msvk.vm_table(dumii-1, dumjj) + (rdumii-dumii) * inv_dum3 *
                  (msvk.vm_table(dumii, dumjj) - msvk.vm_table(dumii-1, dumjj));

                msvk.V_qr(i, k) = (dum1 + (rdumjj - dumjj) * (dum2 - dum1)) * msvk.rhofacr(i, k);

                // number-weighted fall speed:
                dum1 = msvk.vn_table(dumii-1, dumjj-1) + (rdumii-dumii) * inv_dum3 *
                  (msvk.vn_table(dumii, dumjj-1) - msvk.vn_table(dumii-1, dumjj-1));
                dum2 = msvk.vn_table(dumii-1, dumjj) + (rdumii-dumii) * inv_dum3 *
                  (msvk.vn_table(dumii, dumjj) - msvk.vn_table(dumii-1, dumjj));

                msvk.V_nr(i, k) = (dum1 + (rdumjj - dumjj) * (dum2 - dum1)) * msvk.rhofacr(i, k);
              }
              Scalar Co_max_local = msvk.V_qr(i, k) * dt_left * msvk.inv_dzq(i, k);
              if (Co_max_local > lmax) {
                lmax = Co_max_local;
              }
            }, Kokkos::Max<Scalar>(Co_max));

            // compute dt_sub
            int tmpint1 = static_cast<int>(Co_max + 1.0);
            Scalar dt_sub = util::min(dt_left, dt_left / tmpint1);

            int k_temp = (k_qxbot == kbot) ? k_qxbot : (k_qxbot - kdir);

            // calculate fluxes
            util::set_min_max(k_temp, k_qxtop+kdir, kmin, kmax);
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, kmax-kmin+1), [&] (int k_) {
              const int k = kmin + k_;
              msvk.flux_qx(i, k) = msvk.V_qr(i, k) * qr(i, k) * msvk.rho(i, k);
              msvk.flux_nx(i, k) = msvk.V_nr(i, k) * nr(i, k) * msvk.rho(i, k);
            });
            team_member.team_barrier();

            // accumulated precip during time step
            if (k_qxbot == kbot) {
              prt_accum += msvk.flux_qx(i, kbot) * dt_sub;
            }

            Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
              // for top level only (since flux is 0 above)
              int k = k_qxtop;
              // compute flux divergence
              const Scalar fluxdiv_qx = -msvk.flux_qx(i, k) * msvk.inv_dzq(i, k);
              const Scalar fluxdiv_nx = -msvk.flux_nx(i, k) * msvk.inv_dzq(i, k);
              // update prognostic variables
              qr(i, k) += fluxdiv_qx * dt_sub * msvk.inv_rho(i, k);
              nr(i, k) += fluxdiv_nx * dt_sub * msvk.inv_rho(i, k);
            });
            team_member.team_barrier();

            util::set_min_max(k_qxtop - kdir, k_temp, kmin, kmax);
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, kmax-kmin+1), [&] (int k_) {
              const int k = kmin + k_;
              {
                // compute flux divergence
                const Scalar fluxdiv_qx = (msvk.flux_qx(i, k+kdir) - msvk.flux_qx(i, k)) * msvk.inv_dzq(i, k);
                const Scalar fluxdiv_nx = (msvk.flux_nx(i, k+kdir) - msvk.flux_nx(i, k)) * msvk.inv_dzq(i, k);
                // update prognostic variables
                qr(i, k) += fluxdiv_qx * dt_sub * msvk.inv_rho(i, k);
                nr(i, k) += fluxdiv_nx  *dt_sub * msvk.inv_rho(i, k);
              }
            });
            team_member.team_barrier();

            dt_left -= dt_sub;  // update time remaining for sedimentation
            if (k_qxbot != kbot) {
              k_qxbot -= kdir;
            }
          }

          Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
            prt_liq(i) += prt_accum * Globals<Scalar>::INV_RHOW * odt;
          });
        }
      });

    reset(msvk);
  }

};

} // namespace micro_sed
} // namespace p3

#endif
