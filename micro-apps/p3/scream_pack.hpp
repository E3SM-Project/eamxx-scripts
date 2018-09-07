#ifndef INCLUDE_SCREAM_PACK
#define INCLUDE_SCREAM_PACK

namespace scream {
template <int PACKN>
struct Mask {
  enum { packtag = true };
  enum { n = PACKN };

  KOKKOS_FORCEINLINE_FUNCTION Mask (const bool& init) {
    vector_simd for (int i = 0; i < n; ++i) d[i] = init;
  }

  KOKKOS_FORCEINLINE_FUNCTION void set (const int& i, const bool& val) { d[i] = val; }
  KOKKOS_FORCEINLINE_FUNCTION bool operator[] (const int& i) const { return d[i]; }

private:
  char d[n];
};

#define packsimd_gen_assign_op_p(op)                      \
  KOKKOS_FORCEINLINE_FUNCTION                             \
  Pack& operator op (const Pack& a) {                     \
    vector_simd for (int i = 0; i < n; ++i) d[i] op a[i]; \
    return *this;                                         \
  }
#define packsimd_gen_assign_op_s(op)                    \
  KOKKOS_FORCEINLINE_FUNCTION                           \
  Pack& operator op (const scalar& a) {                 \
    vector_simd for (int i = 0; i < n; ++i) d[i] op a;  \
    return *this;                                       \
  }
#define packsimd_gen_assign_op_all(op)          \
  packsimd_gen_assign_op_p(op)                  \
  packsimd_gen_assign_op_s(op)                  \

template <typename SCALAR, int PACKN>
struct Pack {
  enum { packtag = true };
  enum { n = PACKN };

  typedef SCALAR scalar;

  KOKKOS_FORCEINLINE_FUNCTION Pack () {
#ifndef KOKKOS_ENABLE_CUDA
    vector_simd for (int i = 0; i < n; ++i)
      d[i] = std::numeric_limits<scalar>::quiet_NaN();
#endif
  }
  KOKKOS_FORCEINLINE_FUNCTION Pack (const scalar& v) {
    vector_simd for (int i = 0; i < n; ++i) d[i] = v;
  }

  KOKKOS_FORCEINLINE_FUNCTION const scalar& operator[] (const int& i) const { return d[i]; }
  KOKKOS_FORCEINLINE_FUNCTION scalar& operator[] (const int& i) { return d[i]; }

  packsimd_gen_assign_op_all(=)
  packsimd_gen_assign_op_all(+=)
  packsimd_gen_assign_op_all(-=)
  packsimd_gen_assign_op_all(*=)
  packsimd_gen_assign_op_all(/=)

  KOKKOS_FORCEINLINE_FUNCTION void set (const Mask<n>& mask, const scalar& v) {
    vector_simd for (int i = 0; i < n; ++i) if (mask[i]) d[i] = v;
  }
  
private:
  scalar d[n];
};

// Use enable_if and packtag so that we can template on 'Pack' and yet not have
// our operator overloads, in particular, be used for something other than the
// Pack type.
template <typename Pack>
using OnlyPack = typename std::enable_if<Pack::packtag,Pack>::type;
template <typename Pack, typename Return>
using OnlyPackReturn = typename std::enable_if<Pack::packtag,Return>::type;

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> shift_right (const Pack& pm1, const Pack& p) {
  Pack s;
  s[0] = pm1[Pack::n-1];
  vector_simd for (int i = 1; i < Pack::n; ++i) s[i] = p[i-1];
  return s;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> shift_right (const typename Pack::scalar& pm1, const Pack& p) {
  Pack s;
  s[0] = pm1;
  vector_simd for (int i = 1; i < Pack::n; ++i) s[i] = p[i-1];
  return s;
}

#define packsimd_gen_bin_op_pp(op)                                      \
  template <typename Pack> KOKKOS_FORCEINLINE_FUNCTION                  \
  OnlyPack<Pack> operator op (const Pack& a, const Pack& b) {           \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] op b[i];  \
    return c;                                                           \
  }
#define packsimd_gen_bin_op_ps(op)                                      \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyPack<Pack> operator op (const Pack& a, const Scalar& b) {         \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] op b;     \
    return c;                                                           \
  }
#define packsimd_gen_bin_op_sp(op)                                      \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyPack<Pack> operator op (const Scalar& a, const Pack& b) {         \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a op b[i];     \
    return c;                                                           \
  }
#define packsimd_gen_bin_op_all(op)             \
  packsimd_gen_bin_op_pp(op)                    \
  packsimd_gen_bin_op_ps(op)                    \
  packsimd_gen_bin_op_sp(op)

packsimd_gen_bin_op_all(+)
packsimd_gen_bin_op_all(-)
packsimd_gen_bin_op_all(*)
packsimd_gen_bin_op_all(/)

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> cos (const Pack& p) {
  Pack s;
  vector_simd for (int i = 0; i < Pack::n; ++i) s[i] = std::cos(p[i]);
  return s;
}
template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> abs (const Pack& p) {
  Pack s;
  vector_simd for (int i = 0; i < Pack::n; ++i) s[i] = std::abs(p[i]);
  return s;
}
template <typename T, typename pack, typename Scalar> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<pack, Pack<T,pack::n> > min (const Scalar& a, const pack& b) {
  Pack<T,pack::n> c;
  vector_simd for (int i = 0; i < pack::n; ++i) c[i] = ko::min<T>(a, b[i]);
  return c;
}

template <typename T, int n> KOKKOS_INLINE_FUNCTION 
Pack<T,n> pack_range (const T& start) {
  typedef Pack<T,n> pack;
  pack p;
  vector_simd for (int i = 0; i < n; ++i) p[i] = start + i;
  return p;
}

template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack, Mask<Pack::n> > operator >= (const Pack& a, const Scalar& b) {
  Mask<Pack::n> m(false);
  vector_simd for (int i = 0; i < Pack::n; ++i) if (a[i] >= b) m.set(i, true);
  return m;
}
}

#endif
