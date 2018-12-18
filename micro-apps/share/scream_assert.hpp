#ifndef SCREAM_ASSERT_HPP
#define SCREAM_ASSERT_HPP

#include <sstream>
#include <exception>

/*
 * Asserts and error checking for Scream.
 *
 * micro_k* are for error checking within kokkos kernels.
 *
 * Any check with "assert" in the name is disabled for release builds
 *
 * For _msg checks, the msg argument can contain '<<' if not a kernel check.
 */

// Internal do not call directly
#define impl_throw(condition, msg, exception_type)                      \
  do {                                                                  \
    if ( ! (condition) ) {                                              \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition; \
      _ss_ << "\n" << msg;                                              \
      throw exception_type(_ss_.str());                                 \
    }                                                                   \
  } while(0)

#define impl_kthrow(condition, msg)             \
  do {                                          \
    if ( ! (condition) )                        \
      Kokkos::abort(#condition "\n" msg);       \
  } while (0)

#ifndef NDEBUG
#define micro_assert(condition)           impl_throw(condition, "", std::logic_error)
#define micro_kassert(condition)          impl_kthrow(condition, "")
#define micro_assert_msg(condition, msg)  impl_throw(condition, msg, std::logic_error)
#define micro_kassert_msg(condition, msg) impl_kthrow(condition, msg)
#else
#define micro_assert(condition)  ((void) (0))
#define micro_kassert(condition) ((void) (0))
#define micro_assert_msg(condition, msg)  ((void) (0))
#define micro_kassert_msg(condition, msg) ((void) (0))
#endif

#define micro_require(condition)           impl_throw(condition, "", std::logic_error)
#define micro_krequire(condition)          impl_kthrow(condition, "")
#define micro_require_msg(condition, msg)  impl_throw(condition, msg, std::logic_error)
#define micro_krequire_msg(condition, msg) impl_kthrow(condition, msg)

#endif
