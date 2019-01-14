#include "p3_functions_find_impl.hpp"
#include "types.hpp"

namespace p3 {
namespace micro_sed {

/*
 * Explicit instatiation for doing p3 find functions on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace micro_sed
} // namespace p3
