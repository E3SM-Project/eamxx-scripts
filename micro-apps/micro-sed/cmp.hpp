#ifndef INCLUDE_CMP
#define INCLUDE_CMP

#include "types.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

/*
 * A couple micro-app routines for doing comparisons between data arrays. The
 * compare function is the core low-level routine for result comparison between
 * rain-sed implementations.
 */

namespace cmp {

struct TransposeDirection {
  enum Enum { c2f, f2c };
};

// Switch whether i (column index) or k (level index) is the fast
// index. TransposeDirection::c2f makes i faster; f2c makes k faster.
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
static Int compare (const std::string& label, // Label of compared quantities, for use in output
                    // a(1:n) and b(1:n) are compared component-wise.
                    const Scalar* a, const Scalar* b, const Int& n,
                    // The relative error |a(i)-b(i)|/max_i(|a_i|) must be <= tol for success.
                    const Scalar& tol,
                    // If verbose, output to stdout. In any case, return the number of failed components.
                    bool verbose=false) {
  Int nerr = 0;
  Scalar den = 0;
  for (Int i = 0; i < n; ++i)
    den = std::max(den, std::abs(a[i]));
  Scalar worst = 0;
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

static void expect_another_arg (Int i, Int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}

template <typename Scalar>
Scalar compare_property(const std::string& file1_fn, const std::string& file2_fn,
                        const util::FILEPtr& fid1, const util::FILEPtr& fid2, const char* name,
                        bool verbose)
{
  Scalar item1, item2;

  util::read(&item1, 1, fid1);
  util::read(&item2, 1, fid2);

  micro_require_msg(item1 == item2, "Files " << file1_fn << " and " << file2_fn <<
                    " have mismatch in basic property '" << name << "', " << item1 << " != " << item2);

  if (verbose) std::cout << name << " matches, = " << item1 << std::endl;

  return item1;
}

} // namespace cmp

#endif
