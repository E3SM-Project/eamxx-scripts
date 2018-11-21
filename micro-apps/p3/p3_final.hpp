#ifndef MICRO_SED_PACKNOIWS_KOKKOS_HPP
#define MICRO_SED_PACKNOIWS_KOKKOS_HPP

#include "micro_kokkos.hpp"
#include "p3_common.hpp"
#include "scream_pack.hpp"

namespace p3 {
namespace micro_sed {

/*
 * MicroSedFuncFinalKokkos is the implementation of the rain
 * sedimentation component of P3. This implementation is the
 * gold standard for what a Scream class encapsulating a Kokkos
 * kernel(s) should look like.
 *
 * The only point of this class is to encapsulate the micro_sed_func
 * kernel. It is designed to be used as follows:
 *   MicroSedFuncFinalKokkos<Real> msfk(ni, nk);
 *   * initialize inputs (qr, nr, th, dzq, pres) *
 *   * allocate output (prt_liq) *
 *   for (number of time steps)
 *     msfk.micro_sed_func(*args*);
 *
 * In order to support ETI, this class follows the PIMPL pattern.
 */

template <typename Scalar, typename D=DefaultDevice>
struct MicroSedFuncFinalKokkos
{
  //
  // ------- Types --------
  //

  using Pack = BigPack<Scalar>;

  using KT = KokkosTypes<D>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;

  //
  // ------ public API -------
  //

  static constexpr const char* NAME = "final";

  MicroSedFuncFinalKokkos(int num_horz, int num_vert);

  // Rain sed calculation
  //
  // kts: vertical array bound (top)
  // kte: vertical array bound (bottom)
  // its: horizontal array bound
  // ite: horizontal array bound
  // dt: time step
  // qr: rain, mass mixing ratio  (in/out)
  // nr: rain, number mixing ratio (in/out)
  // th: potential temperature                    K
  // dzq: vertical grid spacing                   m
  // pres: pressure                               Pa
  // prt_liq: precipitation rate, total liquid    m s-1  (output)
  void micro_sed_func(
    const Int kts, const Int kte, const int its, const int ite, const Scalar dt,
    const view_2d<Pack>& qr, const view_2d<Pack>& nr,
    const view_2d<const Pack>& th, const view_2d<const Pack>& dzq, const view_2d<const Pack>& pres,
    const view_1d<Scalar>& prt_liq);

  int get_num_vert() const;

  static std::string custom_msg () {
    std::ostringstream out;
    out << " packn=" << SCREAM_PACKN << " small_pack_factor=" << SCREAM_SMALL_PACK_FACTOR;
    return out.str();
  }

  //
  // ---------- Private --------------
  //

#ifndef KOKKOS_ENABLE_CUDA
 private:
#endif

  class Impl;
  std::shared_ptr<Impl> impl;
};

} // namespace micro_sed
} // namespace p3

#ifdef KOKKOS_ENABLE_CUDA
# include "p3_final_impl.hpp"
#endif

#endif
