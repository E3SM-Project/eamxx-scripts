#ifndef INCLUDE_INITIAL_CONDITIONS
#define INCLUDE_INITIAL_CONDITIONS

#include <cassert>
#include <cmath>
#include <vector>

#include "types.hpp"
#include "micro_kokkos.hpp"

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

/* MicroSedData holds data packed by column. */
template <typename Scalar>
struct MicroSedData {
  const Int ni, nk;

protected:
  std::vector<Scalar> buf_;

  void init_ptrs () {
    qr = buf_.data();
    nr = qr + ni*nk;
    th = nr + ni*nk;
    dzq = th + ni*nk;
    pres = dzq + ni*nk;
    prt_liq = pres + ni*nk;
  }

public:
  Real dt;
  bool reverse;
  Real* qr;
  Real* nr;
  Real* th;
  Real* dzq;  // Not sure what 'q' means here, but this is dz [m].
  Real* pres;
  Real* prt_liq;

  MicroSedData(Int ni_, Int nk_)
    : ni(ni_), nk(nk_),
      buf_(5*ni*nk + ni),
      dt(0), reverse(false)
  { init_ptrs(); }

  MicroSedData (const MicroSedData& d)
    : ni(d.ni), nk(d.nk), buf_(d.buf_), dt(d.dt), reverse(d.reverse)
  { init_ptrs(); }

  void copy_data (const MicroSedData& d) {
    assert(ni == d.ni);
    assert(k == d.nk);
    std::copy(d.buf_.begin(), d.buf_.end(), buf_.begin());
  }
};

template <typename Scalar>
void duplicate_columns (MicroSedData<Scalar>& d) {
  const auto copy = [&] (Scalar* v, const Int& i) {
    std::copy(v, v + d.nk, v + i*d.nk);
  };
  for (Int i = 1; i < d.ni; ++i) {
    copy(d.qr, i);
    copy(d.nr, i);
    copy(d.th, i);
    copy(d.dzq, i);
    copy(d.pres, i);
    d.prt_liq[i] = d.prt_liq[0];
  }
}

template <typename Scalar>
void populate (MicroSedData<Scalar>& d, Int kdir) {
  set_hydrostatic(d.nk, d.dzq, d.pres);

  for (Int k = 0; k < d.nk; ++k)
    d.th[k] = consts::th_ref;

  set_rain(d.nk, d.dzq, d.qr, d.nr);

  for (Int i = 0; i < d.ni; ++i)
    d.prt_liq[i] = 0;

  duplicate_columns(d);
  if (kdir == -1) {
    const auto r = reverse_k(d);
    d.copy_data(r);
  }
}

template <typename Scalar>
static MicroSedData<Scalar> reverse_k (const MicroSedData<Scalar>& msd) {
  const auto reverse = [&] (const Scalar* s, Scalar* d) {
    for (Int i = 0; i < msd.ni; ++i) {
      const Scalar* si = s + msd.nk*i;
      Scalar* di = d + msd.nk*i;
      for (Int k = 0; k < msd.nk; ++k)
        di[k] = si[msd.nk - k - 1];
    }
  };
  MicroSedData<Scalar> r(msd);
  r.reverse = ! msd.reverse;
  reverse(msd.qr, r.qr);
  reverse(msd.nr, r.nr);
  reverse(msd.th, r.th);
  reverse(msd.dzq, r.dzq);
  reverse(msd.pres, r.pres);
  return r;
}

} // namespace ic

extern "C" {

void fully_populate_input_data(Int ni, Int nk, Int kdir, Real** qr, Real** nr, Real** th, Real** dzq, Real** pres);

}

#endif
