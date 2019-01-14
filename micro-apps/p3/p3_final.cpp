#include "p3_final_impl.hpp"

namespace p3 {
namespace micro_sed {

/*
 * Explicit instatiation for doing micro-sed (final impl) on Reals using the
 * default device.
 */

template struct MicroSedFuncFinalKokkos<Real, DefaultDevice>;

} // namespace micro_sed
} // namespace p3
