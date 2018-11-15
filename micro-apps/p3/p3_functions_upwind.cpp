#include "p3_functions_upwind.hpp"
#include "types.hpp"

namespace p3 {
namespace micro_sed {

template void Functions<Real,DefaultDevice>::calc_first_order_upwind_step<2>(
  const Unmanaged<view_1d<const Spack> >& rho,
  const Unmanaged<view_1d<const Spack> >& inv_rho,
  const Unmanaged<view_1d<const Spack> >& inv_dzq,
  const MemberType& team,
  const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, const Scalar& dt_sub,
  const view_1d_ptr_array<Spack, 2>& flux,
  const view_1d_ptr_array<Spack, 2>& V,
  const view_1d_ptr_array<Spack, 2>& r);

} // namespace micro_sed
} // namespace p3
