#ifndef P3_CONSTANTS_HPP
#define P3_CONSTANTS_HPP

#include "types.hpp"

namespace p3 {
namespace micro_sed {

/*
 * Mathematical constants used by p3.
 */

template <typename Scalar>
struct Constants
{
  static constexpr Scalar INV_RHOW = 1.e-3;
  static constexpr Scalar RHOW     = 997.0;
  static constexpr Scalar THRD     = 1.0/3.0;
  static constexpr Scalar SXTH     = 1.0/6.0;
  static constexpr Scalar PI       = 3.14159265;
  static constexpr Scalar PIOV6    = PI*SXTH;
  static constexpr Scalar CONS1    = PIOV6*RHOW;
  static constexpr Scalar QSMALL   = 1.e-14;
  static constexpr Scalar NSMALL   = 1.e-16;
  static constexpr Scalar RD       = 287.15;
  static constexpr Scalar RHOSUR   = 100000.0/(RD*273.15);
  static constexpr Scalar CP       = 1005.0;
  static constexpr Scalar INV_CP   = 1.0/CP;
};

template <typename Scalar>
constexpr Scalar Constants<Scalar>::NSMALL;

template <typename Scalar>
struct Globals
{
  static constexpr int VTABLE_DIM0 = 300;
  static constexpr int VTABLE_DIM1 = 10;
  static constexpr int MU_R_TABLE_DIM = 150;

  static vector_2d_t<Scalar> VN_TABLE, VM_TABLE;
  static std::vector<Scalar> MU_R_TABLE;
};

template <typename Scalar>
vector_2d_t<Scalar> Globals<Scalar>::VN_TABLE(VTABLE_DIM0, std::vector<Scalar>(VTABLE_DIM1));

template <typename Scalar>
vector_2d_t<Scalar> Globals<Scalar>::VM_TABLE(VTABLE_DIM0, std::vector<Scalar>(VTABLE_DIM1));

template <typename Scalar>
std::vector<Scalar> Globals<Scalar>::MU_R_TABLE(MU_R_TABLE_DIM);

} // namespace p3
} // namespace micro_sed

#endif
