#ifndef MICRO_SED_VANILLA_HPP
#define MICRO_SED_VANILLA_HPP

#include "array_io.hpp"
#include "util.hpp"
#include "initial_conditions.hpp"

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>

extern "C" {
void p3_init();
Real* c_get_vn_table();
Real* c_get_vm_table();
Real* c_get_mu_r_table();
}

namespace p3 {
namespace micro_sed_vanilla {

// Change this to modify how indices are presented in trace. +1 implies fortran-like indices
#define adjust_indices(index) index + 1

#ifdef TRACE
KOKKOS_INLINE_FUNCTION
void trace_loop(const char* name, int b, int e)
{
  std::printf("%s LOOP %d -> %d\n", name, adjust_indices(b), adjust_indices(e));
}

KOKKOS_INLINE_FUNCTION
void trace_data(const char* name, int i, int k, Real value)
{
  std::printf("%s[%d][%d] = %20.12f\n", name, adjust_indices(i), adjust_indices(k), value);
}
#else
#define trace_loop(n, b, e) (static_cast<void>(0))
#define trace_data(n, i, k, v) (static_cast<void>(0))
#endif

template <typename Real>
struct Globals
{
  static constexpr Real INV_RHOW = 1.e-3;
  static constexpr Real RHOW     = 997.0;
  static constexpr Real THRD     = 1.0/3.0;
  static constexpr Real SXTH     = 1.0/6.0;
  static constexpr Real PI       = 3.14159265;
  static constexpr Real PIOV6    = PI*SXTH;
  static constexpr Real CONS1    = PIOV6*RHOW;
  static constexpr Real QSMALL   = 1.e-14;
  static constexpr Real NSMALL   = 1.e-16;
  static constexpr Real RD       = 287.15;
  static constexpr Real RHOSUR   = 100000.0/(RD*273.15);
  static constexpr Real CP       = 1005.0;
  static constexpr Real INV_CP   = 1.0/CP;

  static vector_2d_t<Real> VN_TABLE, VM_TABLE;
  static std::vector<Real> MU_R_TABLE;
};

template <typename Real>
vector_2d_t<Real> Globals<Real>::VN_TABLE;

template <typename Real>
vector_2d_t<Real> Globals<Real>::VM_TABLE;

template <typename Real>
std::vector<Real> Globals<Real>::MU_R_TABLE;

template <typename Real>
constexpr Real Globals<Real>::NSMALL;

template <typename Real>
void populate_input(const int its, const int ite, const int kts, const int kte,
                    vector_2d_t<Real> & qr, vector_2d_t<Real> & nr, vector_2d_t<Real> & th, vector_2d_t<Real> & dzq, vector_2d_t<Real> & pres, const ic::MicroSedData<Real>* data = nullptr)
{
  const int num_vert = abs(kte - kts) + 1;
  const int num_horz = (ite - its) + 1;

  ic::MicroSedData<Real> default_data(num_horz, num_vert);
  if (data == nullptr) {
    populate(default_data);
    data = &default_data;
  }

  for (int i = 0; i < num_horz; ++i) {
    for (int k = 0; k < num_vert; ++k) {
      qr[i][k]   = data->qr[i*num_vert + k];
      nr[i][k]   = data->nr[i*num_vert + k];
      th[i][k]   = data->th[i*num_vert + k];
      dzq[i][k]  = data->dzq[i*num_vert + k];
      pres[i][k] = data->pres[i*num_vert + k];
    }
  }
}

/**
 * Generate lookup table for rain fallspeed and ventilation parameters
 * the lookup table is two dimensional as a function of number-weighted mean size
 * proportional to qr/Nr and shape parameter mu_r
 */
template <typename Real>
void p3_init_cpp()
{
  static bool is_init = false;
  if (is_init) {
    return;
  }
  is_init = true;

  Globals<Real>::VN_TABLE.resize(300, std::vector<Real>(10));
  Globals<Real>::VM_TABLE.resize(300, std::vector<Real>(10));
  Globals<Real>::MU_R_TABLE.resize(150);

  p3_init();

  Real* vn_table   = c_get_vn_table();
  Real* vm_table   = c_get_vm_table();
  Real* mu_r_table = c_get_mu_r_table();

  for (int i = 0; i < 300; ++i) {
    for (int k = 0; k < 10; ++k) {
      Globals<Real>::VN_TABLE[i][k] = vn_table[300*k + i];
      Globals<Real>::VM_TABLE[i][k] = vm_table[300*k + i];
    }
  }

  for (int i = 0; i < 150; ++i) {
    Globals<Real>::MU_R_TABLE[i] = mu_r_table[i];
  }
}

/**
 * Finds indices in rain lookup table (3)
 */
template <typename Real>
void find_lookupTable_indices_3(int& dumii, int& dumjj, Real& rdumii, Real& rdumjj, Real& inv_dum3,
                                const Real mu_r, const Real lamr)
{
  // find location in scaled mean size space
  Real dum1 = (mu_r+1.) / lamr;
  if (dum1 <= 195.e-6) {
    inv_dum3  = 0.1;
    rdumii = (dum1*1.e6+5.)*inv_dum3;
    rdumii = std::max<Real>(rdumii, 1.);
    rdumii = std::min<Real>(rdumii,20.);
    dumii  = static_cast<int>(rdumii);
    dumii  = std::max(dumii, 1);
    dumii  = std::min(dumii,20);
  }
  else {
    inv_dum3  = Globals<Real>::THRD*0.1;           // i.e. 1/30
    rdumii = (dum1*1.e+6-195.)*inv_dum3 + 20.;
    rdumii = std::max<Real>(rdumii, 20.);
    rdumii = std::min<Real>(rdumii,300.);
    dumii  = static_cast<int>(rdumii);
    dumii  = std::max(dumii, 20);
    dumii  = std::min(dumii,299);
  }

  // find location in mu_r space
  rdumjj = mu_r+1.;
  rdumjj = std::max<Real>(rdumjj,1.);
  rdumjj = std::min<Real>(rdumjj,10.);
  dumjj  = static_cast<int>(rdumjj);
  dumjj  = std::max(dumjj,1);
  dumjj  = std::min(dumjj,9);
}

/**
 * Computes and returns rain size distribution parameters
 */
template <typename Real>
void get_rain_dsd2(const Real qr, Real& nr, Real& mu_r, Real& rdumii, int& dumii, Real& lamr, std::vector<Real> const& mu_r_table,
                   Real& cdistr, Real& logn0r)
{
  if (qr >= Globals<Real>::QSMALL) {
    // use lookup table to get mu
    // mu-lambda relationship is from Cao et al. (2008), eq. (7)

    // find spot in lookup table
    // (scaled N/q for lookup table parameter space_
    nr = std::max(nr, Globals<Real>::NSMALL);
    Real inv_dum = std::pow(qr / (Globals<Real>::CONS1 * nr * 6.0), Globals<Real>::THRD);

    if (inv_dum < 282.e-6) {
      mu_r = 8.282;
    }
    else if (inv_dum >= 282.e-6 && inv_dum < 502.e-6) {
      // interpolate
      rdumii = (inv_dum-250.e-6)*1.e+6*0.5;
      rdumii = std::max<Real>(rdumii,1.0);
      rdumii = std::min<Real>(rdumii,150.0);
      dumii  = static_cast<int>(rdumii);
      dumii  = std::min(149,dumii);
      mu_r   = mu_r_table[dumii-1] + (mu_r_table[dumii] - mu_r_table[dumii-1]) * (rdumii-dumii);
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
void micro_sed_func_vanilla(const int kts, const int kte, const int ni, const int nk, const int its, const int ite, const Real dt,
                            vector_2d_t<Real> & qr, vector_2d_t<Real> & nr,
                            vector_2d_t<Real> const& th, vector_2d_t<Real> const& dzq, vector_2d_t<Real> const& pres,
                            std::vector<Real> & prt_liq)
{
  const int num_vert = abs(kte - kts) + 1;
  const int num_horz = (ite - its) + 1;

  std::vector<Real> V_qr(num_vert), V_nr(num_vert), flux_qx(num_vert), flux_nx(num_vert);

  vector_2d_t<Real> mu_r(num_horz,    std::vector<Real>(num_vert)),
                    lamr(num_horz,    std::vector<Real>(num_vert)),
                    rhofacr(num_horz, std::vector<Real>(num_vert)),
                    inv_dzq(num_horz, std::vector<Real>(num_vert)),
                    rho(num_horz,     std::vector<Real>(num_vert)),
                    inv_rho(num_horz, std::vector<Real>(num_vert)),
                    t(num_horz,       std::vector<Real>(num_vert)),
                    tmparr1(num_horz, std::vector<Real>(num_vert));

  // inverse of thickness of layers
  for (int i = 0; i < num_horz; ++i) {
    for (int k = 0; k < num_vert; ++k) {
      inv_dzq[i][k] = 1 / dzq[i][k];
      t[i][k] = std::pow(pres[i][k] * 1.e-5, Globals<Real>::RD * Globals<Real>::INV_CP) * th[i][k];
    }
  }

  // constants
  const Real odt = 1.0 / dt;

  // direction of vertical leveling
  const int ktop = (kts < kte) ? num_vert-1 : 0;
  const int kbot = (kts < kte) ? 0: num_vert-1;
  const int kdir = (kts < kte) ? 1  : -1;

  // Rain sedimentation:  (adaptivive substepping)
  trace_loop("i_loop_main", 0, num_horz);
  for (int i = 0; i < num_horz; ++i) {

    trace_loop("  k_loop_1", kbot, ktop);
    for (int k = kbot; k != (ktop+kdir); k+=kdir) {
      rho[i][k] = pres[i][k] / (Globals<Real>::RD * t[i][k]);
      inv_rho[i][k] = 1.0 / rho[i][k];
      rhofacr[i][k] = std::pow(Globals<Real>::RHOSUR * inv_rho[i][k], 0.54);
      trace_data("    rhofacr", i, k, rhofacr[i][k]);
    }

    // Note, we are skipping supersaturation checks

    bool log_qxpresent = false;
    int k_qxtop = kbot;

    // find top, determine qxpresent
    for (int k = ktop; k != (kbot-kdir); k-=kdir) {
      if (qr[i][k] >= Globals<Real>::QSMALL) {
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
        if (qr[i][k] >= Globals<Real>::QSMALL) {
          k_qxbot = k;
          break;
        }
      }

      while (dt_left > 1.e-4) {
        Real Co_max = 0.0;
        for (int kk = 0; kk < num_vert; ++kk) {
          V_qr[kk] = 0.0;
          V_nr[kk] = 0.0;
        }

        trace_loop("  k_loop_sedi_r1", k_qxtop, k_qxbot);
        for (int k = k_qxtop; k != (k_qxbot-kdir); k-=kdir) {
          if (qr[i][k] > Globals<Real>::QSMALL) {
            // Compute Vq, Vn:
            nr[i][k] = std::max(nr[i][k], Globals<Real>::NSMALL);
            trace_data("    nr", i, k, nr[i][k]);
            Real rdumii=0.0, tmp1=0.0, tmp2=0.0, rdumjj=0.0, inv_dum3=0.0;
            int dumii=0, dumjj=0;
            get_rain_dsd2(qr[i][k], nr[i][k], mu_r[i][k], rdumii, dumii, lamr[i][k], Globals<Real>::MU_R_TABLE, tmp1, tmp2);
            find_lookupTable_indices_3(dumii, dumjj, rdumii, rdumjj, inv_dum3, mu_r[i][k], lamr[i][k]);

            // mass-weighted fall speed:
            Real dum1 = Globals<Real>::VM_TABLE[dumii-1][dumjj-1] + (rdumii-dumii) * inv_dum3 * \
              (Globals<Real>::VM_TABLE[dumii][dumjj-1] - Globals<Real>::VM_TABLE[dumii-1][dumjj-1]);
            Real dum2 = Globals<Real>::VM_TABLE[dumii-1][dumjj] + (rdumii-dumii) * inv_dum3 * \
              (Globals<Real>::VM_TABLE[dumii][dumjj] - Globals<Real>::VM_TABLE[dumii-1][dumjj]);

            V_qr[k] = (dum1 + (rdumjj - dumjj) * (dum2 - dum1)) * rhofacr[i][k];
            trace_data("    V_qr", 0, k, V_qr[k]);

            // number-weighted fall speed:
            dum1 = Globals<Real>::VN_TABLE[dumii-1][dumjj-1] + (rdumii-dumii) * inv_dum3 * \
              (Globals<Real>::VN_TABLE[dumii][dumjj-1] - Globals<Real>::VN_TABLE[dumii-1][dumjj-1]);
            dum2 = Globals<Real>::VN_TABLE[dumii-1][dumjj] + (rdumii-dumii) * inv_dum3 * \
              (Globals<Real>::VN_TABLE[dumii][dumjj] - Globals<Real>::VN_TABLE[dumii-1][dumjj]);

            V_nr[k] = (dum1 + (rdumjj - dumjj) * (dum2 - dum1)) * rhofacr[i][k];
            trace_data("    V_nr", 0, k, V_nr[k]);
          }
          Co_max = std::max(Co_max, V_qr[k] * dt_left * inv_dzq[i][k]);
          trace_data("  Co_max", 0, 0, Co_max);
        }

        // compute dt_sub
        int tmpint1 = static_cast<int>(Co_max + 1.0);
        Real dt_sub = std::min(dt_left, dt_left / tmpint1);

        int k_temp = (k_qxbot == kbot) ? k_qxbot : (k_qxbot - kdir);

        // calculate fluxes
        trace_loop("  k_flux_loop", k_temp, k_qxtop);
        for (int k = k_temp; k != (k_qxtop+kdir); k+=kdir) {
          flux_qx[k] = V_qr[k] * qr[i][k] * rho[i][k];
          trace_data("    flux_qx", 0, k, flux_qx[k]);
          flux_nx[k] = V_nr[k] * nr[i][k] * rho[i][k];
          trace_data("    flux_nx", 0, k, flux_nx[k]);
        }

        // accumulated precip during time step
        if (k_qxbot == kbot) {
          prt_accum += flux_qx[kbot] * dt_sub;
        }

        // for top level only (since flux is 0 above)
        int k = k_qxtop;
        // compute flux divergence
        Real fluxdiv_qx = -flux_qx[k] * inv_dzq[i][k];
        Real fluxdiv_nx = -flux_nx[k] * inv_dzq[i][k];
        // update prognostic variables
        qr[i][k] += fluxdiv_qx * dt_sub * inv_rho[i][k];
        trace_data("  qr", i, k, qr[i][k]);
        nr[i][k] += fluxdiv_nx * dt_sub * inv_rho[i][k];
        trace_data("  nr", i, k, nr[i][k]);

        trace_loop("  k_flux_div_loop", k_qxtop - kdir, k_temp);
        for (int k = k_qxtop - kdir; k != (k_temp-kdir); k-=kdir) {
          // compute flux divergence
          fluxdiv_qx = (flux_qx[k+kdir] - flux_qx[k]) * inv_dzq[i][k];
          fluxdiv_nx = (flux_nx[k+kdir] - flux_nx[k]) * inv_dzq[i][k];
          // update prognostic variables
          qr[i][k] += fluxdiv_qx * dt_sub * inv_rho[i][k];
          trace_data("    qr", i, k, qr[i][k]);
          nr[i][k] += fluxdiv_nx  *dt_sub * inv_rho[i][k];
          trace_data("    nr", i, k, nr[i][k]);
        }

        dt_left -= dt_sub;  // update time remaining for sedimentation
        if (k_qxbot != kbot) {
          k_qxbot -= kdir;
        }
      }

      trace_data("  prt_liq", i, 0, prt_liq[i]);
      trace_data("  prt_accum", 0, 0, prt_accum);
      prt_liq[i] += prt_accum * Globals<Real>::INV_RHOW * odt;
      trace_data("  prt_liq", i, 0, prt_liq[i]);
    }
  }
}

template <typename Real>
void micro_sed_func_vanilla_wrap(const int kts, const int kte, const int ni, const int nk, const int its, const int ite, const Real dt, const int ts)
{
  const int num_vert = abs(kte - kts) + 1;
  const int num_horz = (ite - its) + 1;

  vector_2d_t<Real> qr(num_horz,    std::vector<Real>(num_vert)),
                    nr(num_horz,    std::vector<Real>(num_vert)),
                    th(num_horz,    std::vector<Real>(num_vert)),
                    dzq(num_horz,   std::vector<Real>(num_vert)),
                    pres(num_horz,  std::vector<Real>(num_vert));

  std::vector<Real> prt_liq(ni);

  std::cout << "Running micro_sed_vanilla with kts=" << kts << ", kte=" << kte << ", ni=" << ni << ", nk=" << nk
            << ", its=" << its << ", ite=" << ite << ", dt=" << dt << ", ts=" << ts << std::endl;

  populate_input(its, ite, kts, kte, qr, nr, th, dzq, pres);

  auto start = std::chrono::steady_clock::now();

  for (int i = 0; i < ts; ++i) {
    micro_sed_func_vanilla<Real>(kts, kte, ni, nk, its, ite, dt, qr, nr, th, dzq, pres, prt_liq);
  }

  auto finish = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

  std::cout << "Time = " << duration.count() / 1000.0 << " seconds." << std::endl;
}

} // namespace p3
} // namespace micro_sed_vanilla

#endif
