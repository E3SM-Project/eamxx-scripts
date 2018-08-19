/* Investigate different vectorization approaches while still getting maximum
   performance on GPU.

   Use 1D advection with a source term to model P3-sedimentation-like code.

   ws: (KH=/ascldap/users/ambradl/lib/kokkos/cpu-dbg; gfortran -g -c -cpp -fopenmp fvec.f90; g++ -std=c++11 vec.cpp fvec.o -I$KH/include -L$KH/lib -lkokkos -ldl -fopenmp)
   blake (SKX): (dp=1; sz=16; KH=/home/ambradl/lib/kokkos/blake; ifort -DVEC_DP=$dp -xcore-avx2 -O3 -fopenmp -cpp -c fvec.f90; icpc -DVEC_PACKN=$sz -DVEC_DP=$dp -DVEC_DEMOTE_M512D -restrict -std=c++11 -xcore-avx2 -O3 vec.cpp fvec.o -I$KH/include -L$KH/lib -lkokkos -ldl -fopenmp -o blake.out)
   bowman (KNL): (dp=1; sz=16; KH=/home/ambradl/lib/kokkos/knl; ifort -cpp -O3 -xmic-avx512 -fopenmp -c -DVEC_DP=$dp fvec.f90; icpc -DVEC_PACKN=$sz -DVEC_DEMOTE_M512D=0 -DVEC_DP=$dp -std=c++11 -restrict -xmic-avx512 -O3 vec.cpp fvec.o -I$KH/include -L$KH/lib -lkokkos -ldl -fopenmp -o bowman.out)
   waterman (P9, V100): (KH=/home/ambradl/lib/waterman-kokkos/gpu; $KH/bin/nvcc_wrapper -DVEC_PACKN=1 -std=c++11 --expt-extended-lambda -fopenmp -O3 vec.cpp -I$KH/include -L$KH/lib -lkokkos -ldl)

   ws lib: (KH=/ascldap/users/ambradl/lib/kokkos/cpu-dbg; g++ -c -std=c++11 vec.cpp -DVEC_LIBRARY -I$KH/include -L$KH/lib -lkokkos -ldl -fopenmp -fPIC; g++ vec.o -shared -fopenmp -o libvec.so)

   To do:
   - workspace ~ #threads, not #cols
 */

#include <Kokkos_Core.hpp>

#include <sys/time.h>

#include <cassert>
#include <cmath>
#include <cstdio>

#ifdef _OPENMP
# include <omp.h>
#endif

#ifndef VEC_FPE
# define VEC_FPE 0
#endif
#if VEC_FPE
# include <xmmintrin.h>
#endif

#ifndef VEC_NCELL
# define VEC_NCELL 128
#endif
#ifndef VEC_PACKN
# define VEC_PACKN 16
#endif
#ifndef VEC_DEMOTE_M512D
# define VEC_DEMOTE_M512D 0
#endif
#ifndef VEC_DP
# define VEC_DP 1
#endif

#if VEC_DP
typedef double Real;
#else
typedef float Real;
#endif
typedef int Int;

#if defined __INTEL_COMPILER
# pragma message "Intel"
# define vector_ivdep _Pragma("ivdep")
# ifdef _OPENMP
#  define vector_simd _Pragma("omp simd")
# else
#  define vector_simd _Pragma("simd")
# endif
#elif defined __GNUG__
# pragma message "GCC"
# define vector_ivdep _Pragma("GCC ivdep")
# define vector_simd _Pragma("GCC ivdep")
# define restrict __restrict__
#else
# define vector_ivdep
# define vector_simd
# define restrict
#endif

extern "C" void f90_step(const int* ncol, const int* nstep, const Real* dt,
                         Real* rho, Real* work);

namespace ko {
#ifdef KOKKOS_ENABLE_CUDA
// Replacements for namespace std functions that don't run on the GPU.

template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }

KOKKOS_INLINE_FUNCTION bool isfinite (const Real& a) {
  return a == a && a != INFINITY && a != -INFINITY;
}

template <typename T> KOKKOS_INLINE_FUNCTION
const T* max_element (const T* const begin, const T* const end) {
  const T* me = begin;
  for (const T* it = begin + 1; it < end; ++it)
    if ( ! (*it < *me)) // use operator<
      me = it;
  return me;
}
#else
using std::min;
using std::max;
using std::isfinite;
using std::max_element;
#endif
} // namespace ko

namespace util {
Real reldif (const Real* a, const Real* b, const Int& n) {
  Real den = 0, num = 0;
  for (Int i = 0; i < n; ++i) {
    den = std::max(den, std::abs(a[i]));
    num = std::max(num, std::abs(a[i] - b[i]));
  }
  return num/den;
}

Real reldif (const Real* a, const Real* b, const Int& n,
             const Int& chunk, const Int& b_ld) {
  Real den = 0, num = 0;
  for (Int i = 0, N = n/chunk; i < N; ++i) {
    for (Int j = 0; j < chunk; ++j) {
      den = std::max(den, std::abs(a[j]));
      num = std::max(num, std::abs(a[j] - b[j]));
    }
    b += b_ld;
  }
  return num/den;
}

static inline double gettime () {
  timeval t;
  gettimeofday(&t, 0);
  static const double us = 1e6;
  return (t.tv_sec*us + t.tv_usec)/us;
}

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

bool eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

void expect_another_arg (Int i, Int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}

// Alloc and dealloc a fixed array. Use this instead of std::vector b/c of
// first-touch issues with std::vector.
template <typename T> class Space {
  T* p_;
public:
  Space (const Int& n) { p_ = new T[n]; }
  ~Space () { delete[] p_; }
  T& operator[] (std::size_t i) { return p_[i]; }
  const T& operator[] (std::size_t i) const { return p_[i]; }
  T* data () const { return p_; }
};

} // namespace util

/* The model problem is 1D constant-in-time, variable-in-space flow with a
   source term. It is intended to be solved using the explicit-Euler,
   first-order upwind method used in P3.
 */
namespace problem {
// nvcc can't handle static const or static constexpr on device, so we need
// #define alternatives.
#define consts_L 2.0
#define consts_xl 0.0
#define consts_xr (consts_xl + L)
#define consts_u_max 0.1
#define consts_ncell VEC_NCELL
#define consts_dx (consts_L/consts_ncell)
#define consts_rho_ref 1.0

struct consts {
  static constexpr Real L     = consts_L;  // m
  static constexpr Real xl    = consts_xl; // m
  static constexpr Real xr    = consts_xr; // m
  static constexpr Real u_max = consts_u_max; // max flow speed, m/s
  static constexpr Int  ncell = consts_ncell; // number of cells
  static constexpr Real dx    = consts_dx;    // width of cell, m
  static constexpr Real rho_ref = consts_rho_ref;
};
constexpr Int consts::ncell;

KOKKOS_INLINE_FUNCTION Real get_x_ctr (const Int& i) {
  return consts_xl + ((i + 0.5)/consts_ncell)*consts_L;
}

template <typename Real>
KOKKOS_INLINE_FUNCTION Real map_x_to_n11 (const Real& x) {
  return 2*(x - consts_xl)/consts_L - 1;
}

// Want something to act like a table.
KOKKOS_INLINE_FUNCTION Real get_src (const Real& x, const Real& rho) {
  const Real rate[] = {0.2, 0.1, 0.05, 0.025, 0.0125};
  constexpr Int tsize = sizeof(rate)/sizeof(*rate) - 1;
  const Real trhodiffmax = 1;
  const auto rhodiff = std::abs(rho - consts::rho_ref);
  if (rhodiff >= trhodiffmax)
    return 0;
  Real src = rhodiff;
  // Throw a min in for safety and to exercise it.
  const auto idx = ko::min<Int>(tsize, tsize*(rhodiff/trhodiffmax));
  return src*rate[idx];
}

template <typename Real>
KOKKOS_INLINE_FUNCTION Real get_u (const Real& x) {
  using std::cos;
  const auto u = consts_u_max*cos(0.2*map_x_to_n11(x) + 0.25);
  //assert(u > 0 && u <= consts::u_max);
  return u;
}

template <typename Real>
KOKKOS_INLINE_FUNCTION Real get_ic (const Real& x) {
  return consts::rho_ref;
}
} // namespace problem

// Reference impl. Straightforward raw arrays with OpenMP across columns.
namespace reference {
// First-order upwind.
void calc_numerical_flux (const Real* rho, Real* flux) {
  using c = problem::consts;
  // Inflow edge.
  flux[0] = c::rho_ref * problem::get_u(c::xl - 0.5*c::dx);
  // Interior and outflow edge.
  for (Int i = 1; i < c::ncell + 1; ++i)
    flux[i] = rho[i-1] * problem::get_u(problem::get_x_ctr(i-1));
}

// work has size ncell+1.
void step (const Real& dt, Real* rho, Real* work) {
  using c = problem::consts;
  auto* flux = work;
  calc_numerical_flux(rho, flux);
  const auto idx = 1/c::dx;
  for (Int i = 0; i < c::ncell; ++i) {
    const auto fluxdiv = idx*(flux[i+1] - flux[i]);
    const auto src = problem::get_src(problem::get_x_ctr(i), rho[i]);
    rho[i] -= dt*(fluxdiv - src);
  }
}

void step (const Int& ncol, const Real& dt, const Int& nstep, Real* rho,
           Real* work) {
  using C = problem::consts;
# pragma omp parallel for
  for (Int c = 0; c < ncol; ++c)
    for (Int i = 0; i < nstep; ++i)
      step(dt, rho + C::ncell*c, work + (C::ncell+1)*c);
}
} // namespace reference

// Assess the cost, if any, of using Kokkos arrays instead of raw arrays.
namespace koarr {
template <typename T>
using Array = Kokkos::View<T*[VEC_NCELL], Kokkos::LayoutRight>;

KOKKOS_INLINE_FUNCTION
void calc_numerical_flux (const Int& c, const Array<Real>& rho, Real& flux_bdy,
                          Array<Real>& flux_interior) {
  using C = problem::consts;
  flux_bdy = C::rho_ref * problem::get_u(C::xl - 0.5*C::dx);
  vector_ivdep for (Int i = 0; i < C::ncell; ++i)
    flux_interior(c,i) = rho(c,i) * problem::get_u(problem::get_x_ctr(i));
}

void step (const Int& c, const Real& dt, Array<Real>& rho,
           Array<Real>& work) {
  using C = problem::consts;
  auto& flux_int = work;
  Real flux_bdy;
  calc_numerical_flux(c, rho, flux_bdy, flux_int);
  const auto idx = 1/C::dx;
  {
    const auto fluxdiv = idx*(flux_int(c,0) - flux_bdy);
    const auto src = problem::get_src(problem::get_x_ctr(0), rho(c,0));
    rho(c,0) -= dt*(fluxdiv - src);
  }
  vector_ivdep for (Int i = 1; i < C::ncell; ++i) {
    const auto fluxdiv = idx*(flux_int(c,i) - flux_int(c,i-1));
    const auto src = problem::get_src(problem::get_x_ctr(i), rho(c,i));
    rho(c,i) -= dt*(fluxdiv - src);
  }
}

void step (const Real& dt, const Int& nstep, Array<Real>& rho,
           Array<Real>& work) {
# pragma omp parallel for
  for (Int c = 0; c < rho.extent_int(0); ++c)  
    for (Int i = 0; i < nstep; ++i)
      step(c, dt, rho, work);
}
} // namespace koarr

// A true Kokkos impl. Works on GPU and uses the threading strategy that seems
// to fit P3 well.
namespace ko {
template <typename T> using Array =
  Kokkos::View<T*[VEC_NCELL], Kokkos::LayoutRight>;
template <typename T> using Array_p1 =
  Kokkos::View<T*[VEC_NCELL+1], Kokkos::LayoutRight>;
template <typename T> using Vector =
  Kokkos::View<T[VEC_NCELL], Kokkos::LayoutRight>;
template <typename T> using Vector_p1 =
  Kokkos::View<T[VEC_NCELL+1], Kokkos::LayoutRight>;

template <typename ExeSpace>
struct OnGpu { enum : bool { value = false }; };
#ifdef KOKKOS_ENABLE_CUDA
template <> struct OnGpu<Kokkos::Cuda> { enum : bool { value = true }; };
#endif

using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
using Member = TeamPolicy::member_type;

KOKKOS_INLINE_FUNCTION
void calc_numerical_flux (const Member& team, const Vector<Real>& rho,
                          const Vector_p1<Real>& flux) {
  Kokkos::single(
    Kokkos::PerTeam(team),
    [&] () {
      flux(0) = consts_rho_ref * problem::get_u(consts_xl - 0.5*consts_dx);
    });
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, consts_ncell),
    [&] (const int& i) {
      flux(i+1) = rho(i) * problem::get_u(problem::get_x_ctr(i));  
    });
  team.team_barrier();
}

KOKKOS_INLINE_FUNCTION
void step (const Member& team, const Real& dt, const Vector<Real>& rho,
           const Vector_p1<Real>& work) {
  auto& flux = work;
  calc_numerical_flux(team, rho, flux);
  const auto idx = 1/consts_dx;
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, consts_ncell),
    [&] (const int& i) {
      const auto fluxdiv = idx*(flux(i+1) - flux(i));
      const auto src = problem::get_src(problem::get_x_ctr(i), rho(i));
      rho(i) -= dt*(fluxdiv - src);
    });
  team.team_barrier();
}

void step (const Real& dt, const Int& nstep, Array<Real>& rho,
           Array_p1<Real>& work) {
  Kokkos::parallel_for(
    TeamPolicy(rho.extent_int(0),
               OnGpu<Kokkos::DefaultExecutionSpace>::value ? VEC_NCELL : 1,
               1),
    KOKKOS_LAMBDA (const Member& team) {
      const auto rhoc = Kokkos::subview(rho, team.league_rank(), Kokkos::ALL);
      const auto workc = Kokkos::subview(work, team.league_rank(), Kokkos::ALL);
      for (Int i = 0; i < nstep; ++i)
        step(team, dt, rhoc, workc);
    });
}
} // namespace ko

// Use pragmas, perhaps helping auto-vectorization. Separate flux on the
// boundary from the interior in case that helps with alignment.
namespace vec1 {
void calc_numerical_flux (const Real* restrict rho, Real& flux_bdy,
                          Real* restrict flux_interior) {
  using c = problem::consts;
  flux_bdy = c::rho_ref * problem::get_u(c::xl - 0.5*c::dx);
  vector_simd for (Int i = 0; i < c::ncell; ++i)
    flux_interior[i] = rho[i] * problem::get_u(problem::get_x_ctr(i));
}

void step (const Real& dt, Real* restrict rho, Real* restrict work) {
  using c = problem::consts;
  auto* flux_int = work;
  Real flux_bdy;
  calc_numerical_flux(rho, flux_bdy, flux_int);
  const auto idx = 1/c::dx;
  {
    const auto neg_flux_div = idx*(flux_bdy - flux_int[0]);
    const auto src = problem::get_src(problem::get_x_ctr(0), rho[0]);
    rho[0] += dt*(src + neg_flux_div);
  }
  vector_ivdep for (Int i = 1; i < c::ncell; ++i) {
    const auto neg_flux_div = idx*(flux_int[i-1] - flux_int[i]);
    const auto src = problem::get_src(problem::get_x_ctr(i), rho[i]);
    rho[i] += dt*(src + neg_flux_div);
  }
}

void step (const Int& ncol, const Real& dt, const Int& nstep, Real* restrict rho,
           Real* restrict work) {
  using C = problem::consts;
# pragma omp parallel for
  for (Int c = 0; c < ncol; ++c)
    for (Int i = 0; i < nstep; ++i)
      step(dt, rho + C::ncell*c, work + C::ncell*c);
}
} // namespace vec1

// Experiment with packs and impls w/in the pack.
namespace vec2 {

#define vec2_gen_assign_op_p(op)                          \
  KOKKOS_FORCEINLINE_FUNCTION                             \
  Pack& operator op (const Pack& a) {                     \
    vector_simd for (int i = 0; i < n; ++i) d[i] op a[i]; \
    return *this;                                         \
  }
#define vec2_gen_assign_op_s(op)                        \
  KOKKOS_FORCEINLINE_FUNCTION                           \
  Pack& operator op (const scalar& a) {                 \
    vector_simd for (int i = 0; i < n; ++i) d[i] op a;  \
    return *this;                                       \
  }
#define vec2_gen_assign_op_all(op)              \
  vec2_gen_assign_op_p(op)                      \
  vec2_gen_assign_op_s(op)                      \

template <typename SCALAR, int PACKN>
struct Pack {
  enum { packtag = true };
  enum { n = PACKN };

  enum { is_pd = std::is_same<SCALAR,double>::value };

  enum { use_avx512f =
#ifdef __AVX512F__
         n % 8 == 0
#else
         false
#endif
  };

  typedef SCALAR scalar;

  KOKKOS_INLINE_FUNCTION Pack () {
#ifndef KOKKOS_ENABLE_CUDA
    vector_simd for (int i = 0; i < n; ++i)
      d[i] = std::numeric_limits<scalar>::quiet_NaN();
#endif
  }
  KOKKOS_INLINE_FUNCTION Pack (const scalar& v) {
    if (use_avx512f && is_pd && n == 8) {
#ifdef __AVX512F__
      as_m512d(0) = _mm512_set_pd(v,v,v,v,v,v,v,v);
#endif
    } else {
      vector_simd for (int i = 0; i < n; ++i) d[i] = v;
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION const scalar& operator[] (const int& i) const { return d[i]; }
  KOKKOS_FORCEINLINE_FUNCTION scalar& operator[] (const int& i) { return d[i]; }

  vec2_gen_assign_op_all(=)
  vec2_gen_assign_op_all(+=)
  vec2_gen_assign_op_all(-=)
  vec2_gen_assign_op_all(*=)
  vec2_gen_assign_op_all(/=)

#ifdef __AVX512F__
  enum { n4 = n/4 };
  enum { n8 = n/8 };
  KOKKOS_FORCEINLINE_FUNCTION const __m512d& as_m512d (const int& i) const
  { return reinterpret_cast<const __m512d*>(d)[i]; }
  KOKKOS_FORCEINLINE_FUNCTION const __m256d& as_m256d (const int& i) const
  { return reinterpret_cast<const __m256d*>(d)[i]; }
  KOKKOS_FORCEINLINE_FUNCTION __m512d& as_m512d (const int& i)
  { return reinterpret_cast<__m512d*>(d)[i]; }
  KOKKOS_FORCEINLINE_FUNCTION __m256d& as_m256d (const int& i)
  { return reinterpret_cast<__m256d*>(d)[i]; }
#endif
  
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

#define vec2_gen_bin_op_pp(op)                                          \
  template <typename Pack> KOKKOS_FORCEINLINE_FUNCTION                  \
  OnlyPack<Pack> operator op (const Pack& a, const Pack& b) {           \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] op b[i];  \
    return c;                                                           \
  }
#define vec2_gen_bin_op_ps(op)                                          \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyPack<Pack> operator op (const Pack& a, const Scalar& b) {         \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] op b;     \
    return c;                                                           \
  }
#define vec2_gen_bin_op_sp(op)                                          \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyPack<Pack> operator op (const Scalar& a, const Pack& b) {         \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a op b[i];     \
    return c;                                                           \
  }
#define vec2_gen_bin_op_all(op)                 \
  vec2_gen_bin_op_pp(op)                        \
  vec2_gen_bin_op_ps(op)                        \
  vec2_gen_bin_op_sp(op)

template <typename Pack> KOKKOS_FORCEINLINE_FUNCTION
OnlyPack<Pack> operator + (const Pack& a, const Pack& b) {
  Pack c;
  if (Pack::use_avx512f && Pack::is_pd && Pack::n == 8) {
#ifdef __AVX512F__
# if VEC_DEMOTE_M512D
    c.as_m256d(0) = _mm256_add_pd(a.as_m256d(0), b.as_m256d(0));
    c.as_m256d(1) = _mm256_add_pd(a.as_m256d(1), b.as_m256d(1));
# else
    c.as_m512d(0) = _mm512_add_pd(a.as_m512d(0), b.as_m512d(0));
# endif
#endif
  } else {
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] + b[i];
  }
  return c;
}
template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION
OnlyPack<Pack> operator + (const Pack& a, const Scalar& b) {
  if (Pack::use_avx512f && Pack::is_pd && Pack::n == 8) {
#ifdef __AVX512F__
    return a + Pack(b);
#endif
  }
  Pack c;
  vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] + b;
  return c;
}
template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION
OnlyPack<Pack> operator + (const Scalar& a, const Pack& b) {
  if (Pack::use_avx512f && Pack::is_pd && Pack::n == 8) {
#ifdef __AVX512F__
    return Pack(a) + b;
#endif
  }
  Pack c;
  vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a + b[i];
  return c;
}

vec2_gen_bin_op_all(-)
vec2_gen_bin_op_all(*)
vec2_gen_bin_op_all(/)

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
  vector_simd for (int i = 0; i < pack::n; ++i) c[i] = std::min<T>(a, b[i]);
  return c;
}

template <typename T, int n>
Pack<T,n> pack_range (const T& start) {
  typedef Pack<T,n> pack;
  pack p;
  if (pack::use_avx512f && pack::is_pd && pack::n == 8) {
#ifdef __AVX512F__
    p.as_m512d(0) = _mm512_set_pd(start+7, start+6, start+5, start+4,
                                  start+3, start+2, start+1, start);
#endif
  } else {
    vector_simd for (int i = 0; i < n; ++i) p[i] = start + i;
  }
  return p;
}

typedef Pack<Real,8> RealPack;
typedef Pack<Int ,8> IntPack;

struct consts : public problem::consts {
  static constexpr Int npack = (ncell + RealPack::n - 1) / RealPack::n;
  static constexpr Int npackp1 = (ncell + RealPack::n) / RealPack::n;
  static constexpr Int packrem = ncell % RealPack::n;
};

RealPack get_x_ctr (const RealPack& i) {
  return consts::xl + ((i + 0.5)/consts::ncell)*consts::L;
}

RealPack get_src (const RealPack& x, const RealPack& rho) {
  using ko::min;
  const Real rate[] = {0.2, 0.1, 0.05, 0.025, 0.0125};
  constexpr Int tsize = sizeof(rate)/sizeof(*rate) - 1;
  const Real trhodiffmax = 1;
  const auto rhodiff = abs(rho - consts::rho_ref);
  RealPack r;
  vector_simd for (Int i = 0; i < RealPack::n; ++i) {
    if (rhodiff[i] >= trhodiffmax) {
      r[i] = 0;
      continue;
    }
    // Throw a min in for safety and to exercise it.
    const auto idx = min<Int>(tsize, tsize*(rhodiff[i]/trhodiffmax));
    r[i] = rhodiff[i]*rate[idx];
  }
  return r;
}

void calc_numerical_flux (const RealPack* restrict rho, Real& flux_bdy,
                          RealPack* restrict flux_interior) {
  using c = consts;
  flux_bdy = c::rho_ref * problem::get_u(c::xl - 0.5*c::dx);
  auto range = pack_range<Real,RealPack::n>(0);
  for (Int i = 0; i < c::npack; ++i) {
    flux_interior[i] = rho[i] * problem::get_u(get_x_ctr(range));
    range += RealPack::n;
  }
}

void step (const RealPack& dt, RealPack* restrict rho, RealPack* restrict work) {
  using c = consts;
  auto* flux_int = work;
  Real flux_bdy;
  calc_numerical_flux(rho, flux_bdy, flux_int);
  const auto idx = 1/c::dx;
  auto range = pack_range<Real,RealPack::n>(0);
  {
    const auto flux_int_im1 = shift_right(flux_bdy, flux_int[0]);
    const auto neg_flux_div = idx*(flux_int_im1 - flux_int[0]);
    const auto src = get_src(get_x_ctr(range), rho[0]);
    rho[0] += dt*(src + neg_flux_div);
  }
  range += RealPack::n;
  for (Int i = 1; i < c::npack; ++i) {
    const auto flux_int_im1 = shift_right(flux_int[i-1], flux_int[i]);
    const auto neg_flux_div = idx*(flux_int_im1 - flux_int[i]);
    const auto src = get_src(get_x_ctr(range), rho[i]);
    rho[i] += dt*(src + neg_flux_div);
    range += RealPack::n;
  }
}

void step (const Int& ncol, const Real& dt, const Int& nstep,
           RealPack* restrict rho, RealPack* restrict work) {
  using C = consts;
# pragma omp parallel for
  for (Int c = 0; c < ncol; ++c)
    for (Int i = 0; i < nstep; ++i)
      step(dt, rho + C::npack*c, work + (C::npack+1)*c);
}
} // namespace vec2

// So far, vec2 shows that intrinsics are not helpful, but packs can be. So do a
// clean impl of packs with no intrinsics.
namespace packsimd {
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

using RealPack = Pack<Real,VEC_PACKN>;
using IntPack = Pack<Int,VEC_PACKN>;

#define consts_npack ((consts_ncell + RealPack::n - 1) / RealPack::n)
struct consts : public problem::consts {
  static constexpr Int npack = consts_npack;
  static constexpr Int npackp1 = (ncell + RealPack::n) / RealPack::n;
  static constexpr Int packrem = ncell % RealPack::n;
};

KOKKOS_INLINE_FUNCTION 
RealPack get_x_ctr (const RealPack& i) {
  return consts_xl + ((i + 0.5)/consts_ncell)*consts_L;
}

KOKKOS_INLINE_FUNCTION 
RealPack get_src (const RealPack& x, const RealPack& rho) {
  using ko::min;
  const Real rate[] = {0.2, 0.1, 0.05, 0.025, 0.0125};
  constexpr Int tsize = sizeof(rate)/sizeof(*rate) - 1;
  const Real trhodiffmax = 1;
  const auto rhodiff = abs(rho - consts_rho_ref);
  RealPack r(0);
  const auto idx = min<Int>(tsize, tsize*(rhodiff/trhodiffmax));
  vector_simd for (Int i = 0; i < RealPack::n; ++i)
    r[i] = rhodiff[i] < trhodiffmax ? rate[idx[i]] : 0;
  r *= rhodiff;
  return r;
}

void calc_numerical_flux (const RealPack* restrict rho, Real& flux_bdy,
                          RealPack* restrict flux_interior) {
  using c = consts;
  flux_bdy = c::rho_ref * problem::get_u(c::xl - 0.5*c::dx);
  auto range = pack_range<Real,RealPack::n>(0);
  for (Int i = 0; i < c::npack; ++i) {
    flux_interior[i] = rho[i] * problem::get_u(get_x_ctr(range));
    range += RealPack::n;
  }
}

void step (const RealPack& dt, RealPack* restrict rho, RealPack* restrict work) {
  using c = consts;
  auto* flux_int = work;
  Real flux_bdy;
  calc_numerical_flux(rho, flux_bdy, flux_int);
  const auto idx = 1/c::dx;
  auto range = pack_range<Real,RealPack::n>(0);
  {
    const auto flux_int_im1 = shift_right(flux_bdy, flux_int[0]);
    const auto neg_flux_div = idx*(flux_int_im1 - flux_int[0]);
    const auto src = get_src(get_x_ctr(range), rho[0]);
    rho[0] += dt*(src + neg_flux_div);
  }
  range += RealPack::n;
  for (Int i = 1; i < c::npack; ++i) {
    const auto flux_int_im1 = shift_right(flux_int[i-1], flux_int[i]);
    const auto neg_flux_div = idx*(flux_int_im1 - flux_int[i]);
    const auto src = get_src(get_x_ctr(range), rho[i]);
    rho[i] += dt*(src + neg_flux_div);
    range += RealPack::n;
  }
}

void step (const Int& ncol, const Real& dt, const Int& nstep,
           RealPack* restrict rho, RealPack* restrict work) {
  using C = consts;
# pragma omp parallel for
  for (Int c = 0; c < ncol; ++c)
    for (Int i = 0; i < nstep; ++i)
      step(dt, rho + C::npack*c, work + C::npackp1*c);
}
} // namespace packsimd

namespace kopack {
using namespace packsimd;

template <typename T> using Array =
  Kokkos::View<Pack<T,RealPack::n>*[consts::npack], Kokkos::LayoutRight>;
template <typename T> using Array_p1 =
  Kokkos::View<Pack<T,RealPack::n>*[consts::npackp1], Kokkos::LayoutRight>;
template <typename T> using Vector =
  Kokkos::View<Pack<T,RealPack::n>[consts::npack], Kokkos::LayoutRight>;
template <typename T> using Vector_p1 =
  Kokkos::View<Pack<T,RealPack::n>[consts::npackp1], Kokkos::LayoutRight>;

template <typename ExeSpace>
struct OnGpu { enum : bool { value = false }; };
#ifdef KOKKOS_ENABLE_CUDA
template <> struct OnGpu<Kokkos::Cuda> { enum : bool { value = true }; };
#endif

using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
using Member = TeamPolicy::member_type;

KOKKOS_INLINE_FUNCTION
void calc_numerical_flux (const Member& team, const Vector<Real>& rho,
                          Real& flux_bdy, const Vector_p1<Real>& flux_interior) {
  Kokkos::single(
    Kokkos::PerTeam(team),
    [&] () {
      flux_bdy = consts_rho_ref * problem::get_u(consts_xl - 0.5*consts_dx);
    });
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, consts_npack),
    [&] (const int& i) {
      const auto range = pack_range<Real,RealPack::n>(RealPack::n*i);
      flux_interior(i) = rho(i) * problem::get_u(get_x_ctr(range));
    });
  team.team_barrier();
}

KOKKOS_INLINE_FUNCTION
void step (const Member& team, const Real& dt, const Vector<Real>& rho,
           const Vector_p1<Real>& work) {
  Real flux_bdy;
  auto& flux_int = work;
  calc_numerical_flux(team, rho, flux_bdy, flux_int);
  const auto idx = 1/consts_dx;
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, consts_npack),
    [&] (const int& i) {
      const auto range = pack_range<Real,RealPack::n>(RealPack::n*i);
      const auto flux_int_im1 =
        i == 0 ?
        shift_right(flux_bdy, flux_int(i)) :
        shift_right(flux_int(i-1), flux_int(i));
      const auto neg_flux_div = idx*(flux_int_im1 - flux_int(i));
      const auto src = get_src(get_x_ctr(range), rho(i));
      rho(i) += dt*(src + neg_flux_div);
    });
  team.team_barrier();
}

void step (const Real& dt, const Int& nstep, Array<Real>& rho,
           Array_p1<Real>& work) {
  Kokkos::parallel_for(
    TeamPolicy(rho.extent_int(0),
               OnGpu<Kokkos::DefaultExecutionSpace>::value ? VEC_NCELL : 1,
               1),
    KOKKOS_LAMBDA (const Member& team) {
      const auto rhoc = Kokkos::subview(rho, team.league_rank(), Kokkos::ALL);
      const auto workc = Kokkos::subview(work, team.league_rank(), Kokkos::ALL);
      for (Int i = 0; i < nstep; ++i)
        step(team, dt, rhoc, workc);
    });
}

namespace cast {
static constexpr int packn = 16, npack = 64;

template <typename T> using Array =
  Kokkos::View<Pack<T,packn>*[npack], Kokkos::LayoutRight>;

// If packn divides ncell, then we can cast to other pack sizes that divide
// packn.
template <typename T, int n> using VectorPack =
  Kokkos::View<Pack<T,n>[npack*(packn/n)], Kokkos::LayoutRight>;

// Or we can case to a scalar.
template <typename T> using VectorScalar =
  Kokkos::View<T[npack*packn], Kokkos::LayoutRight>;

#if 0
template <typename pack, int n> KOKKOS_INLINE_FUNCTION
VectorPack<typename pack::scalar, n> change_pack (pack* v) {
  //static_assert(n <= m && n % m == 0);
  return VectorPack<typename pack::scalar, n>(reinterpret_cast<Pack<typename pack::scalar, n>*>(v));
}
#endif

void cast () {
  Array<Real> base("base", 10);
  const auto bsv = Kokkos::subview(base, 2, Kokkos::ALL);
  VectorPack<Real,4>(reinterpret_cast<Pack<Real,4>*>(bsv.data()));
  //const auto v1 = change_pack<4>(bsv.data()); // equiv to above
  VectorScalar<Real>(reinterpret_cast<Real*>(bsv.data()));
}
}
} // namespace kopack

namespace driver {
Int unittest () {
  using c = problem::consts;
  static constexpr Int ncol = 4, N = c::ncell*ncol;
  static constexpr Real tol = 20*std::numeric_limits<Real>::epsilon();
  const Real dt = 0.5*c::dx/c::u_max;
  const Int nstep = Int(1.2*c::ncell);
  Real ic[N], work[ncol*(c::ncell+1)];
  for (Int c = 0; c < ncol; ++c)
    for (Int i = 0; i < c::ncell; ++i)
      ic[c::ncell*c + i] = problem::get_ic(problem::get_x_ctr(i));
  Real r[N];
  for (Int i = 0; i < N; ++i) r[i] = ic[i];
  reference::step(ncol, dt, nstep, r, work);
  Int nerr = 0;
  Real r1[N];
  {
    for (Int i = 0; i < N; ++i) r1[i] = ic[i];
    f90_step(&ncol, &nstep, &dt, r1, work);
    const auto re = util::reldif(r, r1, N);
    if (re > tol) {
      std::cout << "f90 re " << re << "\n";
      ++nerr;
    }
  }
  {
    for (Int i = 0; i < N; ++i) r1[i] = ic[i];
    vec1::step(ncol, dt, nstep, r1, work);
    if (util::reldif(r, r1, N) > tol) { std::cout << "vec1 failed\n"; ++nerr; }
  }
  {
    vec2::RealPack rho[ncol*vec2::consts::npack],
      work[ncol*vec2::consts::npackp1];
    Real* const rp = reinterpret_cast<Real*>(&rho);
    for (Int i = 0; i < N; ++i) rp[i] = ic[i];
    vec2::step(ncol, dt, nstep, rho, work);
    const auto ld = vec2::consts::npack*vec2::RealPack::n;
    const auto re = util::reldif(r, rp, N, c::ncell, ld);
    if (re > tol) {
      std::cout << "vec2 re " << re << "\n";
      ++nerr;
    }
  }
  {
    packsimd::RealPack rho[ncol*packsimd::consts::npack],
      work[ncol*packsimd::consts::npackp1];
    Real* const rp = reinterpret_cast<Real*>(&rho);
    for (Int i = 0; i < N; ++i) rp[i] = ic[i];
    packsimd::step(ncol, dt, nstep, rho, work);
    const auto ld = packsimd::consts::npack*packsimd::RealPack::n;
    const auto re = util::reldif(r, rp, N, c::ncell, ld);
    if (re > tol) {
      std::cout << "packsimd re " << re << "\n";
      ++nerr;
    }
  }
#ifndef KOKKOS_ENABLE_CUDA
  {
    koarr::Array<Real> rho("rho", ncol), work("work", ncol);
    for (Int c = 0; c < ncol; ++c)
      for (Int i = 0; i < c::ncell; ++i)
        rho(c,i) = ic[c::ncell*c + i];
    koarr::step(dt, nstep, rho, work);
    const auto re = util::reldif(r, &rho(0,0), N);
    if (re > tol) {
      std::cout << "koarr re " << re << "\n";
      ++nerr;
    }
  }
#endif
  {
    ko::Array<Real> rho("rho", ncol);
    ko::Array_p1<Real> work("work", ncol);
    auto rhom = Kokkos::create_mirror_view(rho);
    for (Int c = 0; c < ncol; ++c)
      for (Int i = 0; i < c::ncell; ++i)
        rhom(c,i) = ic[c::ncell*c + i];
    Kokkos::deep_copy(rho, rhom);
    ko::step(dt, nstep, rho, work);
    Kokkos::deep_copy(rhom, rho);
    const auto re = util::reldif(r, &rhom(0,0), N);
    if (re > tol) {
      std::cout << "ko re " << re << "\n";
      ++nerr;
    }
  }
  {
    kopack::Array<Real> rho("rho", ncol);
    kopack::Array_p1<Real> work("work", ncol);
    auto rhom = Kokkos::create_mirror_view(rho);
    Real* const rp = &rhom(0,0)[0];
    for (Int i = 0; i < N; ++i)
      rp[i] = ic[i];
    Kokkos::deep_copy(rho, rhom);
    kopack::step(dt, nstep, rho, work);
    Kokkos::deep_copy(rhom, rho);
    const auto ld = kopack::consts::npack*kopack::RealPack::n;
    const auto re = util::reldif(r, &rhom(0,0)[0], N, c::ncell, ld);
    if (re > tol) {
      std::cout << "kopack re " << re << "\n";
      ++nerr;
    }
  }
  return nerr;
}

void measure_perf (const Int& ncol, const Int& nstep) {
  using C = problem::consts;
  const Int N = C::ncell*ncol;
  const Real dt = 0.5*C::dx/C::u_max;
  double t1, t2;
  util::Space<Real> ic(N), r(N), work(ncol*(C::ncell+1));

  for (Int c = 0; c < ncol; ++c)
    for (Int i = 0; i < C::ncell; ++i)
      ic[C::ncell*c + i] = problem::get_ic(problem::get_x_ctr(i));

#ifndef KOKKOS_ENABLE_CUDA
  for (Int i = 0; i < N; ++i)
    r[i] = ic[i];
  t1 = util::gettime();
  reference::step(ncol, dt, nstep, r.data(), work.data());
  t2 = util::gettime(); printf("ref      %9.3e\n", t2-t1);

  for (Int i = 0; i < N; ++i)
    r[i] = ic[i];
  t1 = util::gettime();
  f90_step(&ncol, &nstep, &dt, r.data(), work.data());
  t2 = util::gettime(); printf("f90      %9.3e\n", t2-t1);

  for (Int i = 0; i < N; ++i)
    r[i] = ic[i];
  t1 = util::gettime();
  vec1::step(ncol, dt, nstep, r.data(), work.data());
  t2 = util::gettime(); printf("vec1     %9.3e\n", t2-t1);

  {
    util::Space<vec2::RealPack> rho(ncol*vec2::consts::npack),
      work(ncol*vec2::consts::npackp1);
    Real* const rp = reinterpret_cast<Real*>(rho.data());
    for (Int i = 0; i < N; ++i) rp[i] = ic[i];
    t1 = util::gettime();
    vec2::step(ncol, dt, nstep, rho.data(), work.data());
    t2 = util::gettime(); printf("vec2     %9.3e\n", t2-t1);
  }
  {
    util::Space<packsimd::RealPack> rho(ncol*packsimd::consts::npack),
      work(ncol*packsimd::consts::npackp1);
    Real* const rp = reinterpret_cast<Real*>(rho.data());
    for (Int i = 0; i < N; ++i) rp[i] = ic[i];
    t1 = util::gettime();
    packsimd::step(ncol, dt, nstep, rho.data(), work.data());
    t2 = util::gettime(); printf("packsimd %9.3e\n", t2-t1);
  }
  {
    koarr::Array<Real> rho("rho", ncol), work("work", ncol);
    for (Int c = 0; c < ncol; ++c)
      for (Int i = 0; i < C::ncell; ++i)
        rho(c,i) = ic[c*C::ncell + i];
    t1 = util::gettime();
    koarr::step(dt, nstep, rho, work);
    t2 = util::gettime(); printf("koarr    %9.3e\n", t2-t1);
  }
#endif
  {
    ko::Array<Real> rho("rho", ncol);
    ko::Array_p1<Real> work("work", ncol);
    auto rhom = Kokkos::create_mirror_view(rho);
    for (Int c = 0; c < ncol; ++c)
      for (Int i = 0; i < C::ncell; ++i)
        rhom(c,i) = ic[c*C::ncell + i];
    Kokkos::deep_copy(rho, rhom);
    t1 = util::gettime();
    ko::step(dt, nstep, rho, work);
    Kokkos::fence();
    t2 = util::gettime(); printf("ko       %9.3e\n", t2-t1);
  }
  {
    kopack::Array<Real> rho("rho", ncol);
    kopack::Array_p1<Real> work("work", ncol);
    auto rhom = Kokkos::create_mirror_view(rho);
    Real* const rp = &rhom(0,0)[0];
    for (Int i = 0; i < N; ++i)
      rp[i] = ic[i];
    Kokkos::deep_copy(rho, rhom);
    t1 = util::gettime();
    kopack::step(dt, nstep, rho, work);
    Kokkos::fence();
    t2 = util::gettime(); printf("kopack   %9.3e\n", t2-t1);
  }
}
} // namespace driver

struct Input {
  Int ncol, nstep;

  Input () : ncol(64), nstep(50000) {}

  void parse (int argc, char** argv) {
    for (Int i = 1; i < argc-1; ++i) {
      if (util::eq(argv[i], "-nc", "--ncol")) {
        util::expect_another_arg(i, argc);
        ncol = std::atoi(argv[++i]);
      } else if (util::eq(argv[i], "-ns", "--nstep")) {
        util::expect_another_arg(i, argc);
        nstep = std::atoi(argv[++i]);
      } else {
        throw std::runtime_error("Unexpected arg.");
      }
    }
  }
};

#ifdef VEC_LIBRARY
extern "C" {
Int vec_get_ncell () { return problem::consts::ncell; }

// img must have ncell*(nstep+1) slots.
void vec_make_image (const Int nstep, Real* const img) {
  using c = problem::consts;

  const Real dt = 0.5*c::dx/c::u_max;

  for (Int i = 0; i < c::ncell; ++i)
    img[i] = problem::get_ic(problem::get_x_ctr(i));

  util::Space<Real> work(c::ncell+1);
  for (Int ti = 0; ti < nstep; ++ti) {
    Real* const prev = img + c::ncell*ti;
    Real* const curr = img + c::ncell*(ti+1);
    std::copy(prev, prev + c::ncell, curr);
    reference::step(1, dt, curr, work.data());
  }
}
} // extern "C"
#else
int main (int argc, char** argv) {
#if VEC_FPE
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
                         ~( _MM_MASK_INVALID |
                            _MM_MASK_DIV_ZERO |
                            _MM_MASK_OVERFLOW |
                            _MM_MASK_UNDERFLOW ));
#endif
  Kokkos::initialize(argc, argv);
  Int nerr = 0;
  do {
    Input in;
    in.parse(argc, argv);
    const Int nthread =
#ifdef _OPENMP
      omp_get_max_threads()
#else
      1
#endif
      ;
    printf("ncell %d pack %d dp %d avx %s demote %d FPE %d nthread %d ncol %d nstep %d\n",
           VEC_NCELL, VEC_PACKN, VEC_DP, util::active_avx_string().c_str(),
           VEC_DEMOTE_M512D, VEC_FPE, nthread, in.ncol, in.nstep);
    nerr = driver::unittest();
    kopack::cast::cast();
    if (nerr) {
      std::cerr << "nerr " << nerr << "\n";
      //break; 
    }
    driver::measure_perf(in.ncol, in.nstep);
  } while (0);
  Kokkos::finalize();
  return nerr;
}
#endif
