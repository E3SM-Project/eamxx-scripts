#ifdef FPE
# include <xmmintrin.h>
#endif

#include "util.hpp"

extern "C" { void dump_arch_f90() { util::dump_arch(); } }

namespace util {

std::string active_avx_string () {
  std::string s;
#if defined __AVX512F__
  s += "-AVX512F";
#endif
#if defined __AVX2__
  s += "-AVX2";
#endif
#if defined __AVX__
  s += "-AVX";
#endif
  return s;
}

void dump_arch () {
  printf("ARCH: dp %d avx %s FPE %d nthread %d\n",
#ifdef DOUBLE_PRECISION
         1,
#else
         0,
#endif
         util::active_avx_string().c_str(),
#ifdef FPE
         1,
#else
         0,
#endif
#ifdef KOKKOS_ENABLE_OPENMP
         Kokkos::OpenMP::concurrency()
#elif defined _OPENMP
         omp_get_max_threads()
#else
         1
#endif
         );
}

void initialize () {
#ifdef FPE
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
                         ~( _MM_MASK_INVALID |
                            _MM_MASK_DIV_ZERO |
                            _MM_MASK_OVERFLOW ));
#endif
}

bool eq (const std::string& a, const char* const b1, const char* const b2) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
    a == std::string("-") + std::string(b1));
}

} // namespace util
