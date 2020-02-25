static const int obj_i = 1;

#ifdef WITH_FPE
#include <xmmintrin.h>

extern "C" void activate_fpe () {
  static unsigned int const exceptions =
    _MM_MASK_INVALID |
    _MM_MASK_DIV_ZERO |
    _MM_MASK_OVERFLOW;
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~exceptions);
}

#else

extern "C" void activate_fpe () {}

#endif
