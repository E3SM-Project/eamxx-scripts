#ifndef INCLUDE_SCREAM_PACK
#define INCLUDE_SCREAM_PACK

//TODO
// - bounds checking define

#include "util.hpp"

namespace scream {
namespace pack {

template <int PACKN>
struct Mask {
  enum { masktag = true };
  enum { n = PACKN };

  KOKKOS_FORCEINLINE_FUNCTION explicit Mask () {}

  KOKKOS_FORCEINLINE_FUNCTION Mask (const bool& init) {
    vector_simd for (int i = 0; i < n; ++i) d[i] = init;
  }

  KOKKOS_FORCEINLINE_FUNCTION void set (const int& i, const bool& val) { d[i] = val; }
  KOKKOS_FORCEINLINE_FUNCTION bool operator[] (const int& i) const { return d[i]; }
  KOKKOS_FORCEINLINE_FUNCTION char& operator[] (const int& i) { return d[i]; }

  bool any () const {
    bool b = false;
    vector_simd for (int i = 0; i < n; ++i) if (d[i]) b = true;
    return b;
  }

private:
  char d[n];
};

template <typename Mask>
using OnlyMask = typename std::enable_if<Mask::masktag,Mask>::type;
template <typename Mask, typename Return>
using OnlyMaskReturn = typename std::enable_if<Mask::masktag,Return>::type;

template <typename Mask> KOKKOS_FORCEINLINE_FUNCTION
OnlyMaskReturn<Mask, void>
loop (const Mask& m, const std::function<void(int)>& f) {
  vector_simd for (int i = 0; i < Mask::n; ++i)
    if (m[i]) f(i);
}

#if 0
#define scream_masked_loop(mask, fn) \
  scream::pack::loop(mask, [&] (const int s) { fn });
#else
#define scream_masked_loop(mask, fn) \
  do { vector_simd for (int s = 0; s < mask.n; ++s) { if (mask[s]) fn } } while (0)
#endif

#define scream_mask_gen_bin_op_mm(op, impl)                   \
  template <typename Mask> KOKKOS_INLINE_FUNCTION             \
  OnlyMask<Mask> operator op (const Mask& a, const Mask& b) { \
    Mask m(false);                                            \
    vector_simd for (int i = 0; i < Mask::n; ++i)             \
      m[i] = a[i] impl b[i];                                  \
    return m;                                                 \
  }

scream_mask_gen_bin_op_mm(&, &&)
scream_mask_gen_bin_op_mm(|, ||)

template <typename Mask> KOKKOS_INLINE_FUNCTION
OnlyMask<Mask> operator ~ (const Mask& m) {
  Mask nm(false);
  vector_simd for (int i = 0; i < Mask::n; ++i) nm[i] = ! m[i];
  return nm;
}

#define scream_pack_gen_assign_op_p(op)                   \
  KOKKOS_FORCEINLINE_FUNCTION                             \
  Pack& operator op (const Pack& a) {                     \
    vector_simd for (int i = 0; i < n; ++i) d[i] op a[i]; \
    return *this;                                         \
  }
#define scream_pack_gen_assign_op_s(op)                 \
  KOKKOS_FORCEINLINE_FUNCTION                           \
  Pack& operator op (const scalar& a) {                 \
    vector_simd for (int i = 0; i < n; ++i) d[i] op a;  \
    return *this;                                       \
  }
#define scream_pack_gen_assign_op_all(op)       \
  scream_pack_gen_assign_op_p(op)               \
  scream_pack_gen_assign_op_s(op)               \

template <typename SCALAR, int PACKN>
struct Pack {
  enum { packtag = true };
  enum { n = PACKN };

  typedef SCALAR scalar;

  KOKKOS_FORCEINLINE_FUNCTION explicit Pack () {
#ifndef KOKKOS_ENABLE_CUDA
    vector_simd for (int i = 0; i < n; ++i)
      d[i] = std::numeric_limits<scalar>::quiet_NaN();
#endif
  }
  KOKKOS_FORCEINLINE_FUNCTION explicit Pack (const scalar& v) {
    vector_simd for (int i = 0; i < n; ++i) d[i] = v;
  }

  template <typename PackIn> KOKKOS_FORCEINLINE_FUNCTION explicit
  Pack (const PackIn& v, typename std::enable_if<PackIn::packtag>::type* = nullptr) {
    static_assert(static_cast<int>(PackIn::n) == static_cast<int>(n),
                  "Pack::n must be the same.");
    vector_simd for (int i = 0; i < n; ++i) d[i] = v[i];
  }

  KOKKOS_FORCEINLINE_FUNCTION const scalar& operator[] (const int& i) const { return d[i]; }
  KOKKOS_FORCEINLINE_FUNCTION scalar& operator[] (const int& i) { return d[i]; }

  scream_pack_gen_assign_op_all(=)
  scream_pack_gen_assign_op_all(+=)
  scream_pack_gen_assign_op_all(-=)
  scream_pack_gen_assign_op_all(*=)
  scream_pack_gen_assign_op_all(/=)

  KOKKOS_FORCEINLINE_FUNCTION
  void set (const Mask<n>& mask, const scalar& v) {
    vector_simd for (int i = 0; i < n; ++i) if (mask[i]) d[i] = v;
  }
  template <typename PackIn> KOKKOS_FORCEINLINE_FUNCTION
  void set (const Mask<n>& mask, const PackIn& p,
            typename std::enable_if<PackIn::packtag>::type* = nullptr) {
    static_assert(static_cast<int>(PackIn::n) == static_cast<int>(n),
                  "Pack::n must be the same.");
    vector_simd for (int i = 0; i < n; ++i) if (mask[i]) d[i] = p[i];
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

// Later, we might support type promotion. For now, caller must explicitly
// promote a pack's scalar type in mixed-type arithmetic.

#define scream_pack_gen_bin_op_pp(op)                                   \
  template <typename Pack> KOKKOS_FORCEINLINE_FUNCTION                  \
  OnlyPack<Pack> operator op (const Pack& a, const Pack& b) {           \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] op b[i];  \
    return c;                                                           \
  }
#define scream_pack_gen_bin_op_ps(op)                                   \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyPack<Pack> operator op (const Pack& a, const Scalar& b) {         \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] op b;     \
    return c;                                                           \
  }
#define scream_pack_gen_bin_op_sp(op)                                   \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyPack<Pack> operator op (const Scalar& a, const Pack& b) {         \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a op b[i];     \
    return c;                                                           \
  }
#define scream_pack_gen_bin_op_all(op)          \
  scream_pack_gen_bin_op_pp(op)                 \
  scream_pack_gen_bin_op_ps(op)                 \
  scream_pack_gen_bin_op_sp(op)

scream_pack_gen_bin_op_all(+)
scream_pack_gen_bin_op_all(-)
scream_pack_gen_bin_op_all(*)
scream_pack_gen_bin_op_all(/)

#define scream_pack_gen_unary_fn(fn, impl)                            \
  template <typename Pack> KOKKOS_INLINE_FUNCTION                     \
  OnlyPack<Pack> fn (const Pack& p) {                                 \
    Pack s;                                                           \
    vector_simd for (int i = 0; i < Pack::n; ++i) s[i] = impl(p[i]);  \
    return s;                                                         \
  }
#define scream_pack_gen_unary_stdfn(fn) scream_pack_gen_unary_fn(fn, std::fn)
scream_pack_gen_unary_stdfn(abs)
scream_pack_gen_unary_stdfn(exp)
scream_pack_gen_unary_stdfn(log)
scream_pack_gen_unary_stdfn(log10)
scream_pack_gen_unary_stdfn(tgamma)

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack, typename Pack::scalar> min (const Pack& p) {
  typename Pack::scalar v(p[0]);
  vector_simd for (int i = 0; i < Pack::n; ++i) v = util::min(v, p[i]);
  return v;
}
template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack, typename Pack::scalar> max (const Pack& p) {
  typename Pack::scalar v(p[0]);
  vector_simd for (int i = 0; i < Pack::n; ++i) v = util::max(v, p[i]);
  return v;
}

#define scream_pack_gen_bin_fn_pp(fn, impl)           \
  template <typename Pack> KOKKOS_INLINE_FUNCTION     \
  OnlyPack<Pack> fn (const Pack& a, const Pack& b) {  \
    Pack s;                                           \
    vector_simd for (int i = 0; i < Pack::n; ++i)     \
      s[i] = impl(a[i], b[i]);                        \
    return s;                                         \
  }
#define scream_pack_gen_bin_fn_ps(fn, impl)                         \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPack<Pack> fn (const Pack& a, const Scalar& b) {              \
    Pack s;                                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      s[i] = impl(a[i], b);                                         \
    return s;                                                       \
  }
#define scream_pack_gen_bin_fn_sp(fn, impl)                         \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPack<Pack> fn (const Scalar& a, const Pack& b) {              \
    Pack s;                                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      s[i] = impl(a, b[i]);                                         \
    return s;                                                       \
  }
#define scream_pack_gen_bin_fn_all(fn, impl)    \
  scream_pack_gen_bin_fn_pp(fn, impl)           \
  scream_pack_gen_bin_fn_ps(fn, impl)           \
  scream_pack_gen_bin_fn_sp(fn, impl)

scream_pack_gen_bin_fn_all(pow, std::pow)
scream_pack_gen_bin_fn_all(min, util::min)
scream_pack_gen_bin_fn_all(max, util::max)

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

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> shift_left (const Pack& pp1, const Pack& p) {
  Pack s;
  s[Pack::n-1] = pp1[0];
  vector_simd for (int i = 0; i < Pack::n-1; ++i) s[i] = p[i+1];
  return s;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> shift_left (const typename Pack::scalar& pp1, const Pack& p) {
  Pack s;
  s[Pack::n-1] = pp1;
  vector_simd for (int i = 0; i < Pack::n-1; ++i) s[i] = p[i+1];
  return s;
}

#define scream_mask_gen_bin_op_pp(op)             \
  template <typename Pack> KOKKOS_INLINE_FUNCTION \
  OnlyPackReturn<Pack, Mask<Pack::n> >            \
  operator op (const Pack& a, const Pack& b) {    \
    Mask<Pack::n> m(false);                       \
    vector_simd for (int i = 0; i < Pack::n; ++i) \
      m[i] = a[i] op b[i];                        \
    return m;                                     \
  }
#define scream_mask_gen_bin_op_ps(op)                               \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPackReturn<Pack, Mask<Pack::n> >                              \
  operator op (const Pack& a, const Scalar& b) {                    \
    Mask<Pack::n> m(false);                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      m[i] = a[i] op b;                                             \
    return m;                                                       \
  }
#define scream_mask_gen_bin_op_sp(op)                               \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPackReturn<Pack, Mask<Pack::n> >                              \
  operator op (const Scalar& a, const Pack& b) {                    \
    Mask<Pack::n> m(false);                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      m[i] = a op b[i];                                             \
    return m;                                                       \
  }
#define scream_mask_gen_bin_op_all(op)          \
  scream_mask_gen_bin_op_pp(op)                 \
  scream_mask_gen_bin_op_ps(op)                 \
  scream_mask_gen_bin_op_sp(op)

scream_mask_gen_bin_op_all(==)
scream_mask_gen_bin_op_all(>=)
scream_mask_gen_bin_op_all(<=)
scream_mask_gen_bin_op_all(>)
scream_mask_gen_bin_op_all(<)

// Index a scalar array with Pack indices, returning a compatible Pack of array
// values.
template<typename Array1, typename IdxPack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<IdxPack, Pack<typename Array1::value_type, IdxPack::n> >
index (const Array1& a, const IdxPack& i0,
       typename std::enable_if<Array1::Rank == 1>::type* = nullptr) {
  Pack<typename Array1::non_const_value_type, IdxPack::n> p;
  vector_simd for (int i = 0; i < IdxPack::n; ++i)
    p[i] = a(i0[i]);
  return p;
}
template<typename Array2, typename IdxPack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<IdxPack, Pack<typename Array2::value_type, IdxPack::n> >
index (const Array2& a, const IdxPack& i0, const IdxPack& i1,
       typename std::enable_if<Array2::Rank == 2>::type* = nullptr) {
  Pack<typename Array2::non_const_value_type, IdxPack::n> p;
  vector_simd for (int i = 0; i < IdxPack::n; ++i)
    p[i] = a(i0[i], i1[i]);
  return p;
}

#if 0
#define scream_maskedpack_gen_assign_op_p(op)                       \
  KOKKOS_FORCEINLINE_FUNCTION                                       \
  Pack& operator op (const Pack& a) {                               \
    vector_simd for (int i = 0; i < n; ++i)                         \
      if(m[i]) d[i] op a[i];                                        \
    return *this;                                                   \
  }
#define scream_maskedpack_gen_assign_op_s(op)   \
  KOKKOS_FORCEINLINE_FUNCTION                   \
  Pack& operator op (const scalar& a) {         \
    vector_simd for (int i = 0; i < n; ++i)     \
      if (m[i]) d[i] op a;                      \
    return *this;                               \
  }
#define scream_maskedpack_gen_assign_op_all(op) \
  scream_maskedpack_gen_assign_op_p(op)         \
  scream_maskedpack_gen_assign_op_s(op)         \

template <typename SCALAR, int PACKN>
struct MaskedPack {
  enum { maskedpacktag = true };
  enum { n = PACKN };

  typedef SCALAR scalar;

  KOKKOS_FORCEINLINE_FUNCTION explicit MaskedPack () {
#ifndef KOKKOS_ENABLE_CUDA
    vector_simd for (int i = 0; i < n; ++i)
      d[i] = std::numeric_limits<scalar>::quiet_NaN();
#endif
  }
  KOKKOS_FORCEINLINE_FUNCTION explicit
  MaskedPack (const Mask<MaskedPack::n>& mask, const scalar& v) {
    vector_simd for (int i = 0; i < n; ++i) m[i] = mask[i];
    vector_simd for (int i = 0; i < n; ++i) if (m[i]) d[i] = v;
  }

  template <typename PackIn> KOKKOS_FORCEINLINE_FUNCTION explicit
  MaskedPack (const PackIn& v,
              typename std::enable_if<PackIn::maskedpacktag>::type* = nullptr) {
    static_assert(static_cast<int>(PackIn::n) == static_cast<int>(n),
                  "MaskPack::n must be the same.");
    vector_simd for (int i = 0; i < n; ++i) d[i] = v.d[i];
    vector_simd for (int i = 0; i < n; ++i) m[i] = v.m[i];
  }

  KOKKOS_FORCEINLINE_FUNCTION const scalar& operator[] (const int& i) const { return d[i]; }
  KOKKOS_FORCEINLINE_FUNCTION scalar& operator[] (const int& i) { return d[i]; }
  KOKKOS_FORCEINLINE_FUNCTION bool mask (const int& i) const { return m[i]; }

  scream_maskedpack_gen_assign_op_all(=)
  scream_maskedpack_gen_assign_op_all(+=)
  scream_maskedpack_gen_assign_op_all(-=)
  scream_maskedpack_gen_assign_op_all(*=)
  scream_maskedpack_gen_assign_op_all(/=)
  
private:
  scalar d[n];
  char m[n];
};

// Use enable_if and packtag so that we can template on 'Pack' and yet not have
// our operator overloads, in particular, be used for something other than the
// Pack type.
template <typename Pack>
using OnlyMaskedPack = typename std::enable_if<Pack::maskedpacktag,Pack>::type;
template <typename Pack, typename Return>
using OnlyMaskedPackReturn = typename std::enable_if<Pack::maskedpacktag,Return>::type;

// Later, we might support type promotion. For now, caller must explicitly
// promote a pack's scalar type in mixed-type arithmetic.

#define scream_maskedpack_gen_bin_op_pp(op)                         \
  template <typename Pack> KOKKOS_FORCEINLINE_FUNCTION              \
  OnlyMaskedPack<Pack> operator op (const Pack& a, const Pack& b) { \
    Pack c;                                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      m[i] = a.mask(i) && a.mask(b);                                \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      if (m[i]) c[i] = a[i] op b[i];                                \
    return c;                                                       \
  }
#define scream_maskedpack_gen_bin_op_ps(op)                             \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyMaskedPack<Pack> operator op (const Pack& a, const Scalar& b) {   \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) m[i] = a.mask(i);     \
    vector_simd for (int i = 0; i < Pack::n; ++i)                       \
      if (m[i]) c[i] = a[i] op b;                                       \
    return c;                                                           \
  }
#define scream_maskedpack_gen_bin_op_sp(op)                             \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyMaskedPack<Pack> operator op (const Scalar& a, const Pack& b) {   \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) m[i] = b.mask(i);     \
    vector_simd for (int i = 0; i < Pack::n; ++i)                       \
      if (m[i]) c[i] = a op b[i];                                       \
    return c;                                                           \
  }
#define scream_maskedpack_gen_bin_op_all(op)    \
  scream_maskedpack_gen_bin_op_pp(op)           \
  scream_maskedpack_gen_bin_op_ps(op)           \
  scream_maskedpack_gen_bin_op_sp(op)

scream_maskedpack_gen_bin_op_all(+)
scream_maskedpack_gen_bin_op_all(-)
scream_maskedpack_gen_bin_op_all(*)
scream_maskedpack_gen_bin_op_all(/)

#define scream_maskedpack_gen_unary_fn(fn, impl)                      \
  template <typename Pack> KOKKOS_INLINE_FUNCTION                     \
  OnlyMaskedPack<Pack> fn (const Pack& p) {                           \
    Pack s;                                                           \
    vector_simd for (int i = 0; i < Pack::n; ++i) s[i] = impl(p[i]);  \
    return s;                                                         \
  }
#define scream_maskedpack_gen_unary_stdfn(fn) scream_maskedpack_gen_unary_fn(fn, std::fn)
scream_maskedpack_gen_unary_stdfn(abs)
scream_maskedpack_gen_unary_stdfn(exp)
scream_maskedpack_gen_unary_stdfn(log)
scream_maskedpack_gen_unary_stdfn(log10)
scream_maskedpack_gen_unary_stdfn(tgamma)

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyMaskedPackReturn<Pack, typename Pack::scalar> min (const Pack& p) {
  typename Pack::scalar v(p[0]);
  vector_simd for (int i = 0; i < Pack::n; ++i) v = util::min(v, p[i]);
  return v;
}
template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyMaskedPackReturn<Pack, typename Pack::scalar> max (const Pack& p) {
  typename Pack::scalar v(p[0]);
  vector_simd for (int i = 0; i < Pack::n; ++i) v = util::max(v, p[i]);
  return v;
}

#define scream_maskedpack_gen_bin_fn_pp(fn, impl)           \
  template <typename Pack> KOKKOS_INLINE_FUNCTION           \
  OnlyMaskedPack<Pack> fn (const Pack& a, const Pack& b) {  \
    Pack s;                                                 \
    vector_simd for (int i = 0; i < Pack::n; ++i)           \
      s[i] = impl(a[i], b[i]);                              \
    return s;                                               \
  }
#define scream_maskedpack_gen_bin_fn_ps(fn, impl)                   \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyMaskedPack<Pack> fn (const Pack& a, const Scalar& b) {        \
    Pack s;                                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      s[i] = impl(a[i], b);                                         \
    return s;                                                       \
  }
#define scream_maskedpack_gen_bin_fn_sp(fn, impl)                   \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyMaskedPack<Pack> fn (const Scalar& a, const Pack& b) {        \
    Pack s;                                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      s[i] = impl(a, b[i]);                                         \
    return s;                                                       \
  }
#define scream_maskedpack_gen_bin_fn_all(fn, impl)  \
  scream_maskedpack_gen_bin_fn_pp(fn, impl)         \
  scream_maskedpack_gen_bin_fn_ps(fn, impl)         \
  scream_maskedpack_gen_bin_fn_sp(fn, impl)

scream_maskedpack_gen_bin_fn_all(pow, std::pow)
scream_maskedpack_gen_bin_fn_all(min, util::min)
scream_maskedpack_gen_bin_fn_all(max, util::max)
#endif

} // namespace pack
} // namespace scream

#endif
