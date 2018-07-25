#ifndef INCLUDE_INITIAL_CONDITIONS
#define INCLUDE_INITIAL_CONDITIONS

#include <cassert>
#include <cmath>

#include "types.hpp"

namespace ic {

struct consts {
  static constexpr Real
    pres_surface = 1e5, // Pa
    pres_top = pres_surface/4, // model top, Pa
    th_ref = 300, // potential temperature reference, K
    g = 9.8, // grav accel, m/s^2
    R_d = 287.058, // specific gas constant for dry air, J/kg/K
    c_p = 1005; // specific heat of air, const pressure, at 300 K, J/kg/K
};

/* z is ordered from surface to top. Determine z mesh by assuming constant
   potential temperature th and hydrostatic pressure p. Combine
     T = th (p/p0)^(R_d/c_p)       T as a function of th
     p = p0 exp(-g z/(R_d T(p)))   hydrostatic pressure profile
   to get (with x = p/p0)
     z = -R_d/g x^(R_d/c_p) log(x) th.
   In the troposphere, we can assume th is roughly constant, so set it to a 80F
   summer day at the surface, 300 K. We are looking at rain, so we don't care
   about anything above ~8km, I think, so set p_top to p_surface/4.
*/
template <typename Scalar>
static void set_hydrostatic (const Int nk, Scalar* const dzq, Scalar* const pres) {
  const auto dp = (consts::pres_top - consts::pres_surface) / nk;
  const auto c1 = consts::R_d/consts::g, c2 = consts::R_d/consts::c_p;
  Real z_prev = 0, z;
  for (Int k = 0; k < nk; ++k) {
    const auto p = consts::pres_surface + (k+1)*dp;
    pres[k] = p - 0.5*dp;
    const auto x = p/consts::pres_surface;
    z = -c1*std::pow(x, c2)*std::log(x)*consts::th_ref;
    assert(z > z_prev); // Monotonically increasing z.
    dzq[k] = z - z_prev;
    assert(k == 0 || dzq[k] > dzq[k-1]); // dz should increase with height.
    z_prev = z;
  }
  assert(z > 8e3); // Make sure we got to over 8km.
}

/* Based on get_rain_dsd2, we want 7.0241e-8 < (qr/nr) < 3.962e-7 to trigger the
   middle condition on inv_dum. So set up a profile that has these values in a
   number of the cells. We'll do a rough unimodal shape whose peak is above this
   range. Thus, all three branches of the condition will be exercised.
*/
template <typename Scalar>
static void set_rain (const Int nk, const Scalar* const dz,
                      Scalar* const qr, Scalar* const nr) {
  static const Real
    nr0 = 1e-2,
    rain_middle = 5e3, // m
    rain_spread = 1e3; // m
  Real z = 0;
  for (Int k = 0; k < nk; ++k) {
    z += dz[k];
    const auto zmid = z - 0.5*dz[k];
    qr[k] = nr[k] = 0;
    if (std::abs(zmid - rain_middle) <= rain_spread) {
      nr[k] = nr0;
      const auto bell = 0.5*(1 + std::cos(2*M_PI*((zmid - rain_middle)/
                                                  (2*rain_spread))));
      qr[k] = nr[k]*bell*1.5*3.962e-7;
    }
  }
}

} // namespace ic

#endif
