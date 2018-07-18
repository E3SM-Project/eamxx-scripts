#ifndef INCLUDE_UTIL
#define INCLUDE_UTIL

#include <cstdio>
#include <sstream>
#include <memory>
#include <exception>

namespace util {

struct FILECloser { void operator() (FILE* fh) { fclose(fh); } };
using FILEPtr = std::unique_ptr<FILE, FILECloser>;

}

#ifndef NDEBUG
#define micro_assert(condition) do {                                    \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
// TODO: Uncomment when Kokkos is in.
// #define micro_kernel_assert(condition) do {     \
//     if ( ! (condition))                         \
//       Kokkos::abort(#condition);                \
//   } while (0)
#else
#define micro_assert(condition)
#define micro_kernel_assert(condition)
#endif
#define micro_throw_if(condition, message) do {                         \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
// TODO: Uncomment when Kokkos is in.
// #define micro_kernel_throw_if(condition, message) do {              \
//     if (condition)                                                  \
//       Kokkos::abort(#condition " led to the exception\n" message);  \
//   } while (0)

#endif
