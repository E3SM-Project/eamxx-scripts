#ifndef MICRO_SED_PACKNOIWS_KOKKOS_HPP
#define MICRO_SED_PACKNOIWS_KOKKOS_HPP

#include "micro_kokkos.hpp"
#include "p3_common.hpp"
#include "scream_pack.hpp"

namespace p3 {
namespace micro_sed {

template <typename Scalar, typename D=DefaultDevice>
struct MicroSedFuncFinalKokkos
{
  //
  // types
  //

  using Pack = BigPack<Scalar>;

  using KT = KokkosTypes<D>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;

  //
  // members
  //

  static constexpr const char* NAME = "final";

private:

  class Impl;
  std::shared_ptr<Impl> impl;

public:
  MicroSedFuncFinalKokkos(int num_horz, int num_vert);

  //TODO This should be made a member function and the first arg dropped.
  static void micro_sed_func(
    MicroSedFuncFinalKokkos& msfk,
    const Int kts, const Int kte, const int its, const int ite, const Scalar dt,
    const view_2d<Pack>& qr, const view_2d<Pack>& nr,
    const view_2d<Pack>& th, const view_2d<Pack>& dzq, const view_2d<Pack>& pres,
    const view_1d<Scalar>& prt_liq);

  int get_num_vert() const;

  static std::string custom_msg () {
    std::ostringstream out;
    out << " packn=" << SCREAM_PACKN << " small_pack_factor=" << SCREAM_SMALL_PACK_FACTOR;
    return out.str();
  }
};

} // namespace micro_sed
} // namespace p3

#ifdef KOKKOS_ENABLE_CUDA
# include "p3_final_impl.hpp"
#endif

#endif
