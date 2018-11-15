#ifndef MICRO_SED_P3_FUNCTIONS_TABLE3_HPP
#define MICRO_SED_P3_FUNCTIONS_TABLE3_HPP

#include "p3_functions.hpp"

namespace p3 {
namespace micro_sed {

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::lookup (const SmallMask<Scalar>& qr_gt_small, Table3& t, const SmallPack<Scalar>& mu_r,
          const SmallPack<Scalar>& lamr) {
  // find location in scaled mean size space
  const auto dum1 = (mu_r+1.) / lamr;
  const auto dum1_lt = qr_gt_small && (dum1 <= 195.e-6);
  if (dum1_lt.any()) {
    scream_masked_loop(dum1_lt, s) {
      const auto inv_dum3 = 0.1;
      auto rdumii = (dum1[s]*1.e6+5.)*inv_dum3;
      rdumii = util::max<Scalar>(rdumii,  1.);
      rdumii = util::min<Scalar>(rdumii, 20.);
      Int dumii = rdumii;
      dumii = util::max(dumii,  1);
      dumii = util::min(dumii, 20);
      t.inv_dum3[s] = inv_dum3;
      t.rdumii[s] = rdumii;
      t.dumii[s] = dumii;
    }
  }
  const auto dum1_gte = qr_gt_small && ! dum1_lt;
  if (dum1_gte.any()) {
    scream_masked_loop(dum1_gte, s) {
      const auto inv_dum3 = Globals<Scalar>::THRD*0.1;
      auto rdumii = (dum1[s]*1.e+6-195.)*inv_dum3 + 20.;
      rdumii = util::max<Scalar>(rdumii, 20.);
      rdumii = util::min<Scalar>(rdumii,300.);
      Int dumii = rdumii;
      dumii = util::max(dumii, 20);
      dumii = util::min(dumii,299);
      t.inv_dum3[s] = inv_dum3;
      t.rdumii[s] = rdumii;
      t.dumii[s] = dumii;
    }
  }

  // find location in mu_r space
  {
    auto rdumjj = mu_r+1.;
    rdumjj = max(rdumjj,1.);
    rdumjj = min(rdumjj,10.);
    IntSmallPack dumjj(rdumjj);
    dumjj  = max(dumjj,1);
    dumjj  = min(dumjj,9);
    t.rdumjj.set(qr_gt_small, rdumjj);
    t.dumjj.set(qr_gt_small, dumjj);
  }
}

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack Functions<S,D>
::apply_table (const SmallMask<Scalar>& qr_gt_small, const view_2d_table& table,
               const Table3& t) {
  const auto rdumii_m_dumii = t.rdumii - Spack(t.dumii);
  const auto t_im1_jm1 = index(table, t.dumii-1, t.dumjj-1);
  const auto dum1 = (t_im1_jm1 + rdumii_m_dumii * t.inv_dum3 *
                     (index(table, t.dumii, t.dumjj-1) - t_im1_jm1));
  const auto t_im1_j = index(table, t.dumii-1, t.dumjj);
  const auto dum2 = (t_im1_j + rdumii_m_dumii * t.inv_dum3 *
                     (index(table, t.dumii, t.dumjj) - t_im1_j));
  return dum1 + (t.rdumjj - Spack(t.dumjj)) * (dum2 - dum1);
}

} // namespace micro_sed
} // namespace p3

#endif
