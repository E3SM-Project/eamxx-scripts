#ifndef MICRO_SED_VANILLA_HPP
#define MICRO_SED_VANILLA_HPP

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>

namespace p3 {
namespace micro_sed_vanilla {

template <typename Real>
using vector_2d_t = std::vector<std::vector<Real> >;

template <typename Real>
struct Consts
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
};

template <typename Real>
void populate_input(const int its, const int ite, const int kts, const int kte,
                    vector_2d_t<Real> & qr, vector_2d_t<Real> & nr, vector_2d_t<Real> & th, vector_2d_t<Real> & dzq, vector_2d_t<Real> & pres)
{
  const int num_vert = kte - kts;
  const int num_horz = ite - its;

  for (int i = 0; i < num_horz; ++i) {
    for (int k = 0; k < num_vert; ++k) {
      qr[i][k]   = 0;
      nr[i][k]   = 0;
      th[i][k]   = 0;
      dzq[i][k]  = 0;
      pres[i][k] = 0;
    }
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
void micro_sed_func_vanilla(const int kts, const int kte, const int ni, const int nk, const int its, const int ite, const int dt,
                            vector_2d_t<Real> & qr, vector_2d_t<Real> & nr,
                            vector_2d_t<Real> const& th, vector_2d_t<Real> const& dzq, vector_2d_t<Real> const& pres,
                            std::vector<Real> & prt_liq)
{
  const int num_vert = kte - kts;
  const int num_horz = ite - its;

  std::vector<Real> V_qr(num_vert), V_nr(num_vert), flux_qx(num_vert), flux_nx(num_vert);

  vector_2d_t<Real> mu_r(num_horz,    std::vector<Real>(num_vert)),
                    lamr(num_horz,    std::vector<Real>(num_vert)),
                    rhofacr(num_horz, std::vector<Real>(num_vert)),
                    inv_dzq(num_horz, std::vector<Real>(num_vert)),
                    rho(num_horz,     std::vector<Real>(num_vert)),
                    inv_rho(num_horz, std::vector<Real>(num_vert)),
                    t(num_horz,       std::vector<Real>(num_vert)),
                    tmparr1(num_horz, std::vector<Real>(num_vert));

  for (int i = 0; i < num_horz; ++i) {
    for (int k = 0; k < num_vert; ++k) {
      inv_dzq[i][k] = 1 / dzq[i][k];
    }
  }

  const Real odt = 1.0 / dt;

  for (int i = 0; i < num_horz; ++i) {
    for (int k = 0; k < num_vert; ++k) {
      tmparr1[i][k] = std::pow(pres[i][k] * 1.e-5, Consts<Real>::RD * Consts<Real>::INV_CP);
    }
  }

}

template <typename Real>
void micro_sed_func_vanilla_wrap(const int kts, const int kte, const int ni, const int nk, const int its, const int ite, const int dt)
{
  const int num_vert = kte - kts;
  const int num_horz = ite - its;

  vector_2d_t<Real> qr(num_horz,    std::vector<Real>(num_vert)),
                    nr(num_horz,    std::vector<Real>(num_vert)),
                    th(num_horz,    std::vector<Real>(num_vert)),
                    dzq(num_horz,   std::vector<Real>(num_vert)),
                    pres(num_horz,  std::vector<Real>(num_vert));

  std::vector<Real> prt_liq(ni);

  populate_input(its, ite, kts, kte, qr, nr, th, dzq, pres);

  auto start = std::chrono::steady_clock::now();

  micro_sed_func_vanilla<Real>(kts, kte, ni, nk, its, ite, dt, qr, nr, th, dzq, pres, prt_liq);

  auto finish = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

  std::cout << "Time = " << duration.count() / 1000.0 << " seconds." << std::endl;
}

} // namespace p3
} // namespace micro_sed_vanilla

#endif
