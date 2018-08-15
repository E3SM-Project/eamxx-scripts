#ifndef MICRO_SED_VANILLA_KOKKOS_HPP
#define MICRO_SED_VANILLA_KOKKOS_HPP

#include "array_io.hpp"
#include "util.hpp"
#include "initial_conditions.hpp"
#include "micro_kokkos.hpp"
#include "micro_sed_vanilla.hpp"

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>

namespace p3 {
namespace micro_sed_vanilla {

template <typename Real>
using kokkos_2d_table_t = Kokkos::View<Real[300][10]>;

template <typename Real>
using kokkos_1d_table_t = Kokkos::View<Real[150]>;

/**
 * Finds indices in rain lookup table (3)
 */
template <typename Real> KOKKOS_FUNCTION
void find_lookupTable_indices_3_kokkos(int& dumii, int& dumjj, Real& rdumii, Real& rdumjj, Real& inv_dum3,
                                       const Real mu_r, const Real lamr)
{
  // find location in scaled mean size space
  Real dum1 = (mu_r+1.) / lamr;
  if (dum1 <= 195.e-6) {
    inv_dum3  = 0.1;
    rdumii = (dum1*1.e6+5.)*inv_dum3;
    rdumii = util::max<Real>(rdumii, 1.);
    rdumii = util::min<Real>(rdumii,20.);
    dumii  = static_cast<int>(rdumii);
    dumii  = util::max(dumii, 1);
    dumii  = util::min(dumii,20);
  }
  else {
    inv_dum3  = Globals<Real>::THRD*0.1;           // i.e. 1/30
    rdumii = (dum1*1.e+6-195.)*inv_dum3 + 20.;
    rdumii = util::max<Real>(rdumii, 20.);
    rdumii = util::min<Real>(rdumii,300.);
    dumii  = static_cast<int>(rdumii);
    dumii  = util::max(dumii, 20);
    dumii  = util::min(dumii,299);
  }

  // find location in mu_r space
  rdumjj = mu_r+1.;
  rdumjj = util::max<Real>(rdumjj,1.);
  rdumjj = util::min<Real>(rdumjj,10.);
  dumjj  = static_cast<int>(rdumjj);
  dumjj  = util::max(dumjj,1);
  dumjj  = util::min(dumjj,9);
}

/**
 * Computes and returns rain size distribution parameters
 */
template <typename Real> KOKKOS_FUNCTION
void get_rain_dsd2_kokkos(const Real qr, Real& nr, Real& mu_r, Real& rdumii, int& dumii, Real& lamr,
                          kokkos_1d_table_t<Real> const& mu_r_table, Real& cdistr, Real& logn0r)
{
  if (qr >= Globals<Real>::QSMALL) {
    // use lookup table to get mu
    // mu-lambda relationship is from Cao et al. (2008), eq. (7)

    // find spot in lookup table
    // (scaled N/q for lookup table parameter space_
    nr = util::max(nr, Globals<Real>::NSMALL);
    Real inv_dum = std::pow(qr / (Globals<Real>::CONS1 * nr * 6.0), Globals<Real>::THRD);

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

    lamr   = std::pow((Globals<Real>::CONS1 *nr *(mu_r+3.0) * (mu_r+2) * (mu_r+1.)/(qr)), Globals<Real>::THRD); // recalculate slope based on mu_r
    Real lammax = (mu_r+1.)*1.e+5;  // check for slope
    Real lammin = (mu_r+1.)*1250.0; // set to small value since breakup is explicitly included (mean size 0.8 mm)

    // apply lambda limiters for rain
    if (lamr < lammin) {
      lamr = lammin;
      nr   = std::exp(3.*std::log(lamr) + std::log(qr) + std::log(std::tgamma(mu_r+1.)) - std::log(std::tgamma(mu_r+4.)))/(Globals<Real>::CONS1);
    }
    else if (lamr > lammax) {
      lamr = lammax;
      nr   = std::exp(3.*std::log(lamr) + std::log(qr) + std::log(std::tgamma(mu_r+1.)) - log(std::tgamma(mu_r+4.)))/(Globals<Real>::CONS1);
    }

    cdistr  = nr/std::tgamma(mu_r+1.);
    logn0r  = std::log10(nr) + (mu_r+1.)*std::log10(lamr) - std::log10(std::tgamma(mu_r+1)); // note: logn0r is calculated as log10(n0r);
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
template <typename Real>
class MicroSedFuncVanillaKokkos
{
 private:
  kokkos_1d_t<Real> V_qr, V_nr, flux_qx, flux_nx;
  kokkos_2d_t<Real> mu_r, lamr, rhofacr, inv_dzq, rho, inv_rho, t, tmparr1;
  kokkos_2d_table_t<Real> vn_table, vm_table;
  kokkos_1d_table_t<Real> mu_r_table;
  int _num_horz, _num_vert;

public:
  MicroSedFuncVanillaKokkos(int num_horz, int num_vert) :
    V_qr("V_qr", num_vert),
    V_nr("V_nr", num_vert),
    flux_qx("flux_qx", num_vert),
    flux_nx("flux_nx", num_vert),
    mu_r("mu_r", num_horz, num_vert),
    lamr("lamr", num_horz, num_vert),
    rhofacr("rhofacr", num_horz, num_vert),
    inv_dzq("inv_dzq", num_horz, num_vert),
    rho("rho", num_horz, num_vert),
    inv_rho("inv_rho", num_horz, num_vert),
    t("t", num_horz, num_vert),
    tmparr1("tmparr1", num_horz, num_vert),
    vn_table("VN_TABLE"), vm_table("VM_TABLE"),
    mu_r_table("MU_R_TABLE"),
    _num_horz(num_horz), _num_vert(num_vert)
  {
    reset();

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

  void reset()
  {
    Kokkos::parallel_for(_num_vert, KOKKOS_LAMBDA(int k) {
      V_qr(k) = 0.0;
      V_nr(k) = 0.0;
      flux_qx(k) = 0.0;
      flux_nx(k) = 0.0;
    });

    Kokkos::parallel_for(_num_horz, KOKKOS_LAMBDA(int i) {
      for (int k = 0; k < _num_vert; ++k) {
        mu_r(i, k)    = 0.0;
        lamr(i, k)    = 0.0;
        rhofacr(i, k) = 0.0;
        inv_dzq(i, k) = 0.0;
        rho(i, k)     = 0.0;
        inv_rho(i, k) = 0.0;
        t(i, k)       = 0.0;
        tmparr1(i, k) = 0.0;
      }
    });
  }

  void micro_sed_func_vanilla_kokkos(const int kts, const int kte, const int ni, const int nk, const int its, const int ite, const Real dt,
                                     kokkos_2d_t<Real> & qr, kokkos_2d_t<Real> & nr,
                                     kokkos_2d_t<Real> const& th, kokkos_2d_t<Real> const& dzq, kokkos_2d_t<Real> const& pres,
                                     kokkos_1d_t<Real> & prt_liq)
  {
    const int num_vert = abs(kte - kts) + 1;
    const int num_horz = (ite - its) + 1;

    // inverse of thickness of layers
    Kokkos::parallel_for(_num_horz, KOKKOS_LAMBDA(int i) {
      for (int k = 0; k < num_vert; ++k) {
        inv_dzq(i, k) = 1 / dzq(i, k);
        t(i, k) = std::pow(pres(i, k) * 1.e-5, Globals<Real>::RD * Globals<Real>::INV_CP) * th(i, k);
      }
    });

    // constants
    const Real odt = 1.0 / dt;

    // direction of vertical leveling
    const int ktop = (kts < kte) ? num_vert-1 : 0;
    const int kbot = (kts < kte) ? 0: num_vert-1;
    const int kdir = (kts < kte) ? 1  : -1;

    // Rain sedimentation:  (adaptivive substepping)
    trace_loop("i_loop_main", 0, num_horz);
    Kokkos::parallel_for(_num_horz, KOKKOS_LAMBDA(int i) {

      trace_loop("  k_loop_1", kbot, ktop);
      for (int k = kbot; k != (ktop+kdir); k+=kdir) {
        rho(i, k) = pres(i, k) / (Globals<Real>::RD * t(i, k));
        inv_rho(i, k) = 1.0 / rho(i, k);
        rhofacr(i, k) = std::pow(Globals<Real>::RHOSUR * inv_rho(i, k), 0.54);
        trace_data("    rhofacr", i, k, rhofacr(i, k));
      }

      // Note, we are skipping supersaturation checks

      bool log_qxpresent = false;
      int k_qxtop = kbot;

      // find top, determine qxpresent
      for (int k = ktop; k != (kbot-kdir); k-=kdir) {
        if (qr(i, k) >= Globals<Real>::QSMALL) {
          log_qxpresent = true;
          k_qxtop = k;
          break;
        }
      }

      // JGF: It appears rain sedimentation is mostly nothing unless log_qxpresent is true
      if (log_qxpresent) {

        Real dt_left = dt;    // time remaining for sedi over full model (mp) time step
        Real prt_accum = 0.0; // precip rate for individual category
        int k_qxbot = 0;

        // find bottom
        for (int k = kbot; k != (k_qxtop+kdir); k+=kdir) {
          if (qr(i, k) >= Globals<Real>::QSMALL) {
            k_qxbot = k;
            break;
          }
        }

        while (dt_left > 1.e-4) {
          Real Co_max = 0.0;
          for (int kk = 0; kk < num_vert; ++kk) {
            V_qr(kk) = 0.0;
            V_nr(kk) = 0.0;
          }

          trace_loop("  k_loop_sedi_r1", k_qxtop, k_qxbot);
          for (int k = k_qxtop; k != (k_qxbot-kdir); k-=kdir) {
            if (qr(i, k) > Globals<Real>::QSMALL) {
              // Compute Vq, Vn:
              nr(i, k) = util::max(nr(i, k), Globals<Real>::NSMALL);
              trace_data("    nr", i, k, nr(i, k));
              Real rdumii=0.0, tmp1=0.0, tmp2=0.0, rdumjj=0.0, inv_dum3=0.0;
              int dumii=0, dumjj=0;
              get_rain_dsd2_kokkos(qr(i, k), nr(i, k), mu_r(i, k), rdumii, dumii, lamr(i, k), mu_r_table, tmp1, tmp2);
              find_lookupTable_indices_3_kokkos(dumii, dumjj, rdumii, rdumjj, inv_dum3, mu_r(i, k), lamr(i, k));

              // mass-weighted fall speed:
              Real dum1 = vm_table(dumii-1, dumjj-1) + (rdumii-dumii) * inv_dum3 * \
                (vm_table(dumii, dumjj-1) - vm_table(dumii-1, dumjj-1));
              Real dum2 = vm_table(dumii-1, dumjj) + (rdumii-dumii) * inv_dum3 * \
                (vm_table(dumii, dumjj) - vm_table(dumii-1, dumjj));

              V_qr(k) = (dum1 + (rdumjj - dumjj) * (dum2 - dum1)) * rhofacr(i, k);
              trace_data("    V_qr", 0, k, V_qr(k));

              // number-weighted fall speed:
              dum1 = vn_table(dumii-1, dumjj-1) + (rdumii-dumii) * inv_dum3 * \
                (vn_table(dumii, dumjj-1) - vn_table(dumii-1, dumjj-1));
              dum2 = vn_table(dumii-1, dumjj) + (rdumii-dumii) * inv_dum3 * \
                (vn_table(dumii, dumjj) - vn_table(dumii-1, dumjj));

              V_nr(k) = (dum1 + (rdumjj - dumjj) * (dum2 - dum1)) * rhofacr(i, k);
              trace_data("    V_nr", 0, k, V_nr(k));
            }
            Co_max = util::max(Co_max, V_qr(k) * dt_left * inv_dzq(i, k));
            trace_data("  Co_max", 0, 0, Co_max);
          }

          // compute dt_sub
          int tmpint1 = static_cast<int>(Co_max + 1.0);
          Real dt_sub = util::min(dt_left, dt_left / tmpint1);

          int k_temp = (k_qxbot == kbot) ? k_qxbot : (k_qxbot - kdir);

          // calculate fluxes
          trace_loop("  k_flux_loop", k_temp, k_qxtop);
          for (int k = k_temp; k != (k_qxtop+kdir); k+=kdir) {
            flux_qx(k) = V_qr(k) * qr(i, k) * rho(i, k);
            trace_data("    flux_qx", 0, k, flux_qx(k));
            flux_nx(k) = V_nr(k) * nr(i, k) * rho(i, k);
            trace_data("    flux_nx", 0, k, flux_nx(k));
          }

          // accumulated precip during time step
          if (k_qxbot == kbot) {
            prt_accum += flux_qx(kbot) * dt_sub;
          }

          // for top level only (since flux is 0 above)
          int k = k_qxtop;
          // compute flux divergence
          Real fluxdiv_qx = -flux_qx(k) * inv_dzq(i, k);
          Real fluxdiv_nx = -flux_nx(k) * inv_dzq(i, k);
          // update prognostic variables
          qr(i, k) += fluxdiv_qx * dt_sub * inv_rho(i, k);
          trace_data("  qr", i, k, qr(i, k));
          nr(i, k) += fluxdiv_nx * dt_sub * inv_rho(i, k);
          trace_data("  nr", i, k, nr(i, k));

          trace_loop("  k_flux_div_loop", k_qxtop - kdir, k_temp);
          for (int k = k_qxtop - kdir; k != (k_temp-kdir); k-=kdir) {
            // compute flux divergence
            fluxdiv_qx = (flux_qx(k+kdir) - flux_qx(k)) * inv_dzq(i, k);
            fluxdiv_nx = (flux_nx(k+kdir) - flux_nx(k)) * inv_dzq(i, k);
            // update prognostic variables
            qr(i, k) += fluxdiv_qx * dt_sub * inv_rho(i, k);
            trace_data("    qr", i, k, qr(i, k));
            nr(i, k) += fluxdiv_nx  *dt_sub * inv_rho(i, k);
            trace_data("    nr", i, k, nr(i, k));
          }

          dt_left -= dt_sub;  // update time remaining for sedimentation
          if (k_qxbot != kbot) {
            k_qxbot -= kdir;
          }
        }

        trace_data("  prt_liq", i, 0, prt_liq(i));
        trace_data("  prt_accum", 0, 0, prt_accum);
        prt_liq(i) += prt_accum * Globals<Real>::INV_RHOW * odt;
        trace_data("  prt_liq", i, 0, prt_liq(i));
      }
    });

    reset();
  }
};

template <typename Real>
void populate_kokkos_from_vec(const int num_horz, const int num_vert, vector_2d_t<Real> const& vec, kokkos_2d_t<Real>& device)
{
  typename kokkos_2d_t<Real>::HostMirror mirror = Kokkos::create_mirror_view(device);

  for (int i = 0; i < num_horz; ++i) {
    for (int k = 0; k < num_vert; ++k) {
      mirror(i, k) = vec[i][k];
    }
  }

  Kokkos::deep_copy(device, mirror);

}

template <typename Real>
void micro_sed_func_vanilla_kokkos_wrap(const int kts, const int kte, const int ni, const int nk, const int its, const int ite, const Real dt, const int ts)
{
  const int num_vert = abs(kte - kts) + 1;
  const int num_horz = (ite - its) + 1;

  vector_2d_t<Real> qr_v(num_horz,    std::vector<Real>(num_vert)),
                    nr_v(num_horz,    std::vector<Real>(num_vert)),
                    th_v(num_horz,    std::vector<Real>(num_vert)),
                    dzq_v(num_horz,   std::vector<Real>(num_vert)),
                    pres_v(num_horz,  std::vector<Real>(num_vert));

  std::cout << "Running micro_sed_vanilla_kokkos with kts=" << kts << ", kte=" << kte << ", ni=" << ni << ", nk=" << nk
            << ", its=" << its << ", ite=" << ite << ", dt=" << dt << ", ts=" << ts << std::endl;

  populate_input(its, ite, kts, kte, qr_v, nr_v, th_v, dzq_v, pres_v);

  kokkos_2d_t<Real> qr("qr", num_horz, num_vert),
    nr("nr", num_horz, num_vert),
    th("th", num_horz, num_vert),
    dzq("dzq", num_horz, num_vert),
    pres("pres", num_horz, num_vert);

  kokkos_1d_t<Real> prt_liq("prt_liq", ni);

  for (auto item : { std::make_pair(&qr_v, &qr), std::make_pair(&nr_v, &nr), std::make_pair(&th_v, &th),
        std::make_pair(&dzq_v, &dzq), std::make_pair(&pres_v, &pres)}) {
    populate_kokkos_from_vec(num_horz, num_vert, *(item.first), *(item.second));
  }

  Kokkos::parallel_for(ni, KOKKOS_LAMBDA(int i) {
    prt_liq(i) = 0.0;
  });

  MicroSedFuncVanillaKokkos<Real> msvk(num_horz, num_vert);

  auto start = std::chrono::steady_clock::now();

  for (int i = 0; i < ts; ++i) {
    msvk.micro_sed_func_vanilla_kokkos(kts, kte, ni, nk, its, ite, dt, qr, nr, th, dzq, pres, prt_liq);
  }

  auto finish = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

  std::cout << "Time = " << duration.count() / 1000.0 << " seconds." << std::endl;
}

} // namespace p3
} // namespace micro_sed_vanilla

#endif
