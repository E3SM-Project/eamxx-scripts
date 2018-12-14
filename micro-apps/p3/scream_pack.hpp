#ifndef INCLUDE_SCREAM_PACK
#define INCLUDE_SCREAM_PACK

//TODO
// - bounds checking define

#include "util.hpp"
#include "scream_arch.hpp"

namespace scream {
namespace pack {

/*
 * API for using "packed" data in scream. Packs are just bundles of N
 * scalars wthin a single object. Using packed data makes it much easier
 * to get good vectorization with C++.
 */

template <int PACKN>
struct Mask {
  // One tends to think a short boolean type would be useful here, but that is
  // bad for vectorization. int or long are best.
  typedef long type;

  enum { masktag = true };
  enum { n = PACKN };

  KOKKOS_FORCEINLINE_FUNCTION explicit Mask () {}

  KOKKOS_FORCEINLINE_FUNCTION Mask (const bool& init) {
    //vector_simd // Intel 18 is having an issue with this loop.
    for (int i = 0; i < n; ++i) d[i] = init;
  }

  KOKKOS_FORCEINLINE_FUNCTION void set (const int& i, const bool& val) { d[i] = val; }
  KOKKOS_FORCEINLINE_FUNCTION bool operator[] (const int& i) const { return d[i]; }

  KOKKOS_FORCEINLINE_FUNCTION bool any () const {
    bool b = false;
    vector_simd for (int i = 0; i < n; ++i) if (d[i]) b = true;
    return b;
  }

  KOKKOS_FORCEINLINE_FUNCTION bool all () const {
    bool b = true;
    vector_simd for (int i = 0; i < n; ++i) if ( ! d[i]) b = false;
    return b;
  }

private:
  type d[n];
};

template <typename Mask>
using OnlyMask = typename std::enable_if<Mask::masktag,Mask>::type;
template <typename Mask, typename Return>
using OnlyMaskReturn = typename std::enable_if<Mask::masktag,Return>::type;

#define scream_masked_loop(mask, s)                         \
  vector_simd for (int s = 0; s < mask.n; ++s) if (mask[s])

#define scream_mask_gen_bin_op_mm(op, impl)                   \
  template <typename Mask> KOKKOS_INLINE_FUNCTION             \
  OnlyMask<Mask> operator op (const Mask& a, const Mask& b) { \
    Mask m(false);                                            \
    vector_simd for (int i = 0; i < Mask::n; ++i)             \
      if (a[i] impl b[i]) m.set(i, true);                     \
    return m;                                                 \
  }

scream_mask_gen_bin_op_mm(&&, &&)
scream_mask_gen_bin_op_mm(||, ||)

template <typename Mask> KOKKOS_INLINE_FUNCTION
OnlyMask<Mask> operator ! (const Mask& m) {
  Mask nm(false);
  vector_simd for (int i = 0; i < Mask::n; ++i) nm.set(i, ! m[i]);
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

  typedef typename std::remove_const<SCALAR>::type scalar;

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

  template <typename PackIn> KOKKOS_FORCEINLINE_FUNCTION explicit
  Pack (const Mask<Pack::n>& m, const PackIn& p) {
    static_assert(static_cast<int>(PackIn::n) == static_cast<int>(n),
                  "Pack::n must be the same.");
#ifndef KOKKOS_ENABLE_CUDA
    vector_simd for (int i = 0; i < n; ++i)
      d[i] = m[i] ? p[i] : std::numeric_limits<scalar>::quiet_NaN();
#else
    vector_simd for (int i = 0; i < n; ++i)
      d[i] = m[i] ? p[i] : 0;
#endif
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

// max(init, max(p(mask)))
template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack, typename Pack::scalar>
max (const Mask<Pack::n>& mask, typename Pack::scalar init, const Pack& p) {
  vector_simd for (int i = 0; i < Pack::n; ++i)
    if (mask[i]) init = util::max(init, p[i]);
  return init;
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
      s[i] = impl<typename Pack::scalar>(a[i], b);                  \
    return s;                                                       \
  }
#define scream_pack_gen_bin_fn_sp(fn, impl)                         \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPack<Pack> fn (const Scalar& a, const Pack& b) {              \
    Pack s;                                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      s[i] = impl<typename Pack::scalar>(a, b[i]);                  \
    return s;                                                       \
  }
#define scream_pack_gen_bin_fn_all(fn, impl)    \
  scream_pack_gen_bin_fn_pp(fn, impl)           \
  scream_pack_gen_bin_fn_ps(fn, impl)           \
  scream_pack_gen_bin_fn_sp(fn, impl)

scream_pack_gen_bin_fn_all(min, util::min)
scream_pack_gen_bin_fn_all(max, util::max)

// On Intel 17 for KNL, I'm getting a ~1-ulp diff on const Scalar& b. I don't
// understand its source. But, in any case, I'm writing a separate impl here to
// get around that.
//scream_pack_gen_bin_fn_all(pow, std::pow)
template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> pow (const Pack& a, const Scalar/*&*/ b) {
  Pack s;
  vector_simd for (int i = 0; i < Pack::n; ++i)
    s[i] = std::pow<typename Pack::scalar>(a[i], b);
  return s;
}

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
      if (a[i] op b[i]) m.set(i, true);           \
    return m;                                     \
  }
#define scream_mask_gen_bin_op_ps(op)                               \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPackReturn<Pack, Mask<Pack::n> >                              \
  operator op (const Pack& a, const Scalar& b) {                    \
    Mask<Pack::n> m(false);                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      if (a[i] op b) m.set(i, true);                                \
    return m;                                                       \
  }
#define scream_mask_gen_bin_op_sp(op)                               \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPackReturn<Pack, Mask<Pack::n> >                              \
  operator op (const Scalar& a, const Pack& b) {                    \
    Mask<Pack::n> m(false);                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      if (a op b[i]) m.set(i, true);                                \
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
OnlyPackReturn<IdxPack, Pack<typename Array1::non_const_value_type, IdxPack::n> >
index (const Array1& a, const IdxPack& i0,
       typename std::enable_if<Array1::Rank == 1>::type* = nullptr) {
  Pack<typename Array1::non_const_value_type, IdxPack::n> p;
  vector_simd for (int i = 0; i < IdxPack::n; ++i)
    p[i] = a(i0[i]);
  return p;
}
template<typename Array2, typename IdxPack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<IdxPack, Pack<typename Array2::non_const_value_type, IdxPack::n> >
index (const Array2& a, const IdxPack& i0, const IdxPack& i1,
       typename std::enable_if<Array2::Rank == 2>::type* = nullptr) {
  Pack<typename Array2::non_const_value_type, IdxPack::n> p;
  vector_simd for (int i = 0; i < IdxPack::n; ++i)
    p[i] = a(i0[i], i1[i]);
  return p;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack,Int> npack(const Int& nscalar) {
  return (nscalar + Pack::n - 1) / Pack::n;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> range (const typename Pack::scalar& start) {
  typedef Pack pack;
  pack p;
  vector_simd for (int i = 0; i < Pack::n; ++i) p[i] = start + i;
  return p;
}

template <typename T>
using BigPack = Pack<T, SCREAM_PACKN>;
template <typename T>
using SmallPack = Pack<T, SCREAM_PACKN / SCREAM_SMALL_PACK_FACTOR>;

using RealPack = BigPack<Real>;
using IntPack = BigPack<Int>;
using RealSmallPack = SmallPack<Real>;
using IntSmallPack = SmallPack<Int>;

template <typename T, typename ...Parms, int pack_size> KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<Kokkos::View<T**, Parms...> >
scalarize (const Kokkos::View<Pack<T, pack_size>**, Parms...>& vp) {
  return Unmanaged<Kokkos::View<T**, Parms...> >(
    reinterpret_cast<T*>(vp.data()), vp.extent_int(0), pack_size * vp.extent_int(1));
}

template <typename T, typename ...Parms, int pack_size> KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<Kokkos::View<T*, Parms...> >
scalarize (const Kokkos::View<Pack<T, pack_size>*, Parms...>& vp) {
  return Unmanaged<Kokkos::View<T*, Parms...> >(
    reinterpret_cast<T*>(vp.data()), pack_size * vp.extent_int(0));
}

template <int new_pack_size,
          typename T, typename ...Parms, int old_pack_size>
KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<Kokkos::View<Pack<T, new_pack_size>**, Parms...> >
repack (const Kokkos::View<Pack<T, old_pack_size>**, Parms...>& vp) {
  static_assert(new_pack_size > 0 &&
                old_pack_size % new_pack_size == 0,
                "New pack size must divide old pack size.");
  return Unmanaged<Kokkos::View<Pack<T, new_pack_size>**, Parms...> >(
    reinterpret_cast<Pack<T, new_pack_size>*>(vp.data()),
    vp.extent_int(0),
    (old_pack_size / new_pack_size) * vp.extent_int(1));
}

// shrinking
template <int new_pack_size,
          typename T, typename ...Parms, int old_pack_size,
          typename std::enable_if<(old_pack_size >= new_pack_size), int>::type = 0>
KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<Kokkos::View<Pack<T, new_pack_size>*, Parms...> >
repack (const Kokkos::View<Pack<T, old_pack_size>*, Parms...>& vp) {
  static_assert(new_pack_size > 0 &&
                old_pack_size % new_pack_size == 0,
                "New pack size must divide old pack size.");
  return Unmanaged<Kokkos::View<Pack<T, new_pack_size>*, Parms...> >(
    reinterpret_cast<Pack<T, new_pack_size>*>(vp.data()),
    (old_pack_size / new_pack_size) * vp.extent_int(0));
}

// growing
template <int new_pack_size,
          typename T, typename ...Parms, int old_pack_size,
          typename std::enable_if<(old_pack_size < new_pack_size), int>::type = 0>
KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<Kokkos::View<Pack<T, new_pack_size>*, Parms...> >
repack (const Kokkos::View<Pack<T, old_pack_size>*, Parms...>& vp) {
  static_assert(new_pack_size > 0 &&
                new_pack_size % old_pack_size == 0,
                "Old pack size must divide new pack size.");
  micro_kassert(vp.extent_int(0) % (new_pack_size / old_pack_size) == 0);
  return Unmanaged<Kokkos::View<Pack<T, new_pack_size>*, Parms...> >(
    reinterpret_cast<Pack<T, new_pack_size>*>(vp.data()),
    vp.extent_int(0) / (new_pack_size / old_pack_size));
}

template <typename T, typename ...Parms> KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<Kokkos::View<SmallPack<T>**, Parms...> >
smallize (const Kokkos::View<BigPack<T>**, Parms...>& vp) {
  return repack<SCREAM_PACKN / SCREAM_SMALL_PACK_FACTOR>(vp);
}

template <typename T, typename ...Parms> KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<Kokkos::View<SmallPack<T>*, Parms...> >
smallize (const Kokkos::View<BigPack<T>*, Parms...>& vp) {
  return repack<SCREAM_PACKN / SCREAM_SMALL_PACK_FACTOR>(vp);
}

} // namespace pack
} // namespace scream

#endif
