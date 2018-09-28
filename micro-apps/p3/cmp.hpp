#ifndef INCLUDE_CMP
#define INCLUDE_CMP

#include "types.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

namespace cmp {

struct TransposeDirection {
  enum Enum { c2f, f2c };
};

template <TransposeDirection::Enum direction, typename Scalar>
void transpose(const Scalar* sv, Scalar* dv, Int ni, Int nk) {
  for (Int k = 0; k < nk; ++k)
    for (Int i = 0; i < ni; ++i)
      if (direction == TransposeDirection::c2f)
        dv[ni*k + i] = sv[nk*i + k];
      else
        dv[nk*i + k] = sv[ni*k + i];
};

template <typename Scalar>
static Int compare (const std::string& label, const Scalar* a,
                    const Scalar* b, const Int& n, const Real& tol, bool verbose=false) {
  Int nerr = 0;
  Real den = 0;
  for (Int i = 0; i < n; ++i)
    den = std::max(den, std::abs(a[i]));
  Real worst = 0;
  for (Int i = 0; i < n; ++i) {
    const auto num = std::abs(a[i] - b[i]);
    if (num > tol*den ||
        std::isnan(a[i]) || std::isinf(a[i]) ||
        std::isnan(b[i]) || std::isinf(b[i])) {
      ++nerr;
      if (verbose) {
        std::cout << label << " bad idx: " << i << std::fixed << std::setprecision(12)
                  << std::setw(20) << a[i] << " " << b[i] << std::endl;
      }
      worst = std::max(worst, num);
    }
  }
  if (nerr)
    std::cout << label << " nerr " << nerr << " worst " << (worst/den)
              << " with denominator " << den << "\n";
  return nerr;
}

} // namespace cmp

#endif
