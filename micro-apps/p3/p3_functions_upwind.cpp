#include "p3_functions_upwind_impl.hpp"
#include "types.hpp"

namespace p3 {
namespace micro_sed {

#define ETI_UPWIND(nfield)                                              \
  template void Functions<Real,DefaultDevice>                           \
  ::calc_first_order_upwind_step<nfield>(                               \
    const Unmanaged<view_1d<const Spack> >& rho,                        \
    const Unmanaged<view_1d<const Spack> >& inv_rho,                    \
    const Unmanaged<view_1d<const Spack> >& inv_dzq,                    \
    const MemberType& team,                                             \
    const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, \
    const Scalar& dt_sub,                                               \
    const view_1d_ptr_array<Spack, nfield>& flux,                       \
    const view_1d_ptr_array<Spack, nfield>& V,                          \
    const view_1d_ptr_array<Spack, nfield>& r);
ETI_UPWIND(2)
#undef ETI_UPWIND

} // namespace micro_sed
} // namespace p3
