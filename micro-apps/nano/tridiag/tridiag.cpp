#include <sys/time.h>
#include <cassert>
#include <limits>

#include <Kokkos_Core.hpp>

#ifdef TRIDIAG_HAVE_CUSPARSE
# include <cusparse.h>
#endif
#ifdef TRIDIAG_HAVE_MKL
# include <mkl.h>
#endif

/*
  A = [0 b0 c0; a1 b1 c1; ...; an bn 0]
  X fast index is j

  todo:
  x get test1 working on gpu
  x get a perf test going
  x try shared mem
  x impl pcr
  x make a perf test that covers thomas, cr, pcr
  x try LayoutStride, esp with pcr. nope: LayoutStride has overhead
    (revealed by order = {2,1,0}), so it's unusable.
  x make a cusparse perf test
  x 2x2 opt suggests i should try thomas in cr, after all. no.
  x cr_a1x1: hand-code dl,d,u layout since LayoutStride won't cut it.
  x switch interface to take dl,d,du
  x homme gpu version that packs A differently. also can put A and X updates
    together. with this version, nrhs = 1 case should call homme version. then
    check perf against cusparse nrhs=1 case.
  x raw cuda
  x perf test for cr_amxm vs cusparse
  x then opt cr_amxm, following a1x?? examples
  x call thomas with pack
  x make a blas perf test: (gttrf, gttrs); (dttrfb, dttrsb) if available
  x opt thomas
  - make an analysis routine that returns a small struct of POD indicating which
    alg and with what parms to run for a given set of prob parms. then make a
    top-level routine that switches based on the struct.
 */

#if defined __INTEL_COMPILER
# define vector_ivdep _Pragma("ivdep")
# ifdef _OPENMP
#  define vector_simd _Pragma("omp simd")
# else
#  define vector_simd _Pragma("simd")
# endif
#elif defined __GNUG__
# define vector_ivdep _Pragma("GCC ivdep")
# define vector_simd _Pragma("GCC ivdep")
# define restrict __restrict__
#else
# define vector_ivdep
# define vector_simd
# define restrict
#endif

namespace util {
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

void expect_another_arg (int i, int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}
} // namespace util

namespace ko {
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }
} // namespace ko

// Just the Pack goodies we need.
namespace pack {
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
    // Quiet NaNs don't work on Cuda.
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

private:
  scalar d[n];
};

template <typename Pack>
using OnlyPack = typename std::enable_if<Pack::packtag,Pack>::type;
template <typename Pack, typename Return>
using OnlyPackReturn = typename std::enable_if<Pack::packtag,Return>::type;

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
} // namespace pack

template <typename TeamMember>
KOKKOS_INLINE_FUNCTION
int get_thread_id_within_team (const TeamMember& team) {
  return team.team_rank();
}

template <typename TeamMember>
KOKKOS_INLINE_FUNCTION
int get_team_nthr (const TeamMember& team) {
  return team.team_size();
}

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
Kokkos::TeamPolicy<ExeSpace>
get_default_team_policy (const int nprob, const int nrow, const int nrhs) {
#ifdef KOKKOS_ENABLE_OPENMP
# ifdef TRIDIAG_MIMIC_GPU
  return Kokkos::TeamPolicy<ExeSpace>(nprob, omp_get_max_threads(), 1);
# else
  const int per = (omp_get_max_threads() + nprob - 1)/nprob;
  return Kokkos::TeamPolicy<ExeSpace>(nprob, per, 1);
# endif
#else
  return Kokkos::TeamPolicy<ExeSpace>(nprob, 1, 1);
#endif
}

#ifdef KOKKOS_ENABLE_CUDA
KOKKOS_INLINE_FUNCTION
int get_thread_id_within_team (const Kokkos::Impl::CudaTeamMember& team) {
#ifdef __CUDA_ARCH__
  // Can't use team.team_rank() here because vector direction also uses physical
  // threads but TeamMember types don't expose that information.
  return blockDim.x * threadIdx.y + threadIdx.x;
#else
  assert(0);
  return -1;
#endif
}

KOKKOS_INLINE_FUNCTION
int get_team_nthr (const Kokkos::Impl::CudaTeamMember& team) {
#ifdef __CUDA_ARCH__
  return blockDim.x * blockDim.y;
#else
  assert(0);
  return -1;
#endif
}

// This is not intended to be optimal for the tridiag problem; rather,
// it is intended to reflect what a column physics param would likely
// use.
template <>
Kokkos::TeamPolicy<Kokkos::Cuda>
get_default_team_policy<Kokkos::Cuda> (const int nprob, const int nrow, const int nrhs) {
  const int nwarp_per_gpu = 2560;
  const int tpw = 32;
  const int nwarp_per_team =
    std::min(16, // max warp/team
             std::max(1,
                      std::max((nrow + tpw - 1)/tpw,
                               (nwarp_per_gpu + nprob - 1)/nprob)));
  // switch to (nprob, nwarp_per_gpu, tpw) to test vec
  return Kokkos::TeamPolicy<Kokkos::Cuda>(nprob, tpw*nwarp_per_team, 1);
}
#endif

// NB: The caller must provide the team_barrier after this function
// returns before A is accessed.
template <typename TeamMember, typename TridiagDiag>
KOKKOS_INLINE_FUNCTION
void thomas_factorize (const TeamMember& team,
                       TridiagDiag dl, TridiagDiag d, TridiagDiag du) {
  const int nrow = d.extent_int(0);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  const auto f = [&] () {
    for (int i = 1; i < nrow; ++i) {
      dl(i) /= d(i-1);
      d (i) -= dl(i) * du(i-1);
    }
  };
  Kokkos::single(Kokkos::PerTeam(team), f);
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas_solve (const TeamMember& team,
                   TridiagDiag dl, TridiagDiag d, TridiagDiag du,
                   DataArray X) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  const auto f = [&] (const int& j) {
    const auto g = [&] () {
      for (int i = 1; i < nrow; ++i)
        X(i,j) -= dl(i) * X(i-1,j);
      X(nrow-1,j) /= d(nrow-1);
      for (int i = nrow-1; i > 0; --i)
        X(i-1,j) = (X(i-1,j) - du(i-1) * X(i,j)) / d(i-1);
    };
    Kokkos::single(Kokkos::PerThread(team), g);
  };
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nrhs), f);
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas (const TeamMember& team,
             TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X) {
  const int nrow = d.extent_int(0);
  assert(X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  thomas_factorize(team, dl, d, du);
  team.team_barrier();
  thomas_solve(team, dl, d, du, X);
}

template <typename DT, typename XT>
KOKKOS_INLINE_FUNCTION
void thomas_a1x1 (DT* const dl, DT* d, DT* const du, XT* X, const int nrow) {
  for (int i = 1; i < nrow; ++i) {
    const auto dli = dl[i] / d[i-1];
    d[i] -= dli * du[i-1];
    X[i] -= dli * X[i-1];
  }
  X[nrow-1] /= d[nrow-1];
  for (int i = nrow-1; i > 0; --i)
    X[i-1] = (X[i-1] - du[i-1] * X[i]) / d[i-1];
}

template <typename DT, typename XT>
KOKKOS_INLINE_FUNCTION
void thomas_a1xm (DT* const dl, DT* d, DT* const du, XT* X,
                  const int nrow, const int nrhs) {
  for (int i = 1; i < nrow; ++i) {
    const auto dli = dl[i] / d[i-1];
    d[i] -= dli * du[i-1];
    auto* const xim1 = X + (i-1)*nrhs;
    auto* const xi = X + i*nrhs;
    for (int j = 0; j < nrhs; ++j)
      xi[j] -= dli * xim1[j];
  }
  {
    auto* const xi = X + (nrow-1)*nrhs;
    for (int j = 0; j < nrhs; ++j)
      xi[j] /= d[nrow-1];
  }
  for (int i = nrow-1; i > 0; --i) {
    auto* const xim1 = X + (i-1)*nrhs;
    auto* const xi = X + i*nrhs;
    for (int j = 0; j < nrhs; ++j)
      xim1[j] = (xim1[j] - du[i-1] * xi[j]) / d[i-1];
  }
}

template <typename DT, typename XT>
KOKKOS_INLINE_FUNCTION
void thomas_amxm (DT* const dl, DT* d, DT* const du, XT* X,
                  const int nrow, const int nrhs) {
  for (int i = 1; i < nrow; ++i) {
    const int ios = i*nrhs;
    const int im1os = (i-1)*nrhs;
    auto* const dli = dl + ios;
    auto* const di = d + ios;
    auto* const dim1 = d + im1os;
    auto* const duim1 = du + im1os;
    auto* const xim1 = X + im1os;
    auto* const xi = X + ios;
    for (int j = 0; j < nrhs; ++j) {
      const auto dlij = dli[j] / dim1[j];
      di[j] -= dlij * duim1[j];
      xi[j] -= dlij * xim1[j];
    }
  }
  {
    const int ios = (nrow-1)*nrhs;
    auto* const di = d + ios;
    auto* const xi = X + ios;
    for (int j = 0; j < nrhs; ++j)
      xi[j] /= di[j];
  }
  for (int i = nrow-1; i > 0; --i) {
    const int ios = i*nrhs;
    const int im1os = (i-1)*nrhs;
    auto* const dim1 = d + im1os;
    auto* const duim1 = du + im1os;
    auto* const xim1 = X + im1os;
    auto* const xi = X + ios;
    for (int j = 0; j < nrhs; ++j)
      xim1[j] = (xim1[j] - duim1[j] * xi[j]) / dim1[j];
  }
}

template <typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas (TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
             typename std::enable_if<TridiagDiag::rank == 1>::type* = 0,
             typename std::enable_if<DataArray::rank == 1>::type* = 0,
             typename std::enable_if<std::is_same<typename DataArray::array_layout,
                                                  Kokkos::LayoutRight>::value ||
                                     std::is_same<typename DataArray::array_layout,
                                                  Kokkos::LayoutLeft>::value>::type* = 0) {
  const int nrow = d.extent_int(0);
  assert( X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  thomas_a1x1(dl.data(), d.data(), du.data(), X.data(), nrow);
}

template <typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas (TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
             typename std::enable_if<TridiagDiag::rank == 1>::type* = 0,
             typename std::enable_if<DataArray::rank == 2>::type* = 0,
             typename std::enable_if<std::is_same<typename DataArray::array_layout,
                                                  Kokkos::LayoutRight>::value>::type* = 0) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert( X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  thomas_a1xm(dl.data(), d.data(), du.data(), X.data(), nrow, nrhs);
}

template <typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas (TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
             typename std::enable_if<TridiagDiag::rank == 2>::type* = 0,
             typename std::enable_if<DataArray::rank == 2>::type* = 0) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(X .extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(dl.extent_int(1) == nrhs);
  assert(d .extent_int(1) == nrhs);
  assert(du.extent_int(1) == nrhs);
  thomas_amxm(dl.data(), d.data(), du.data(), X.data(), nrow, nrhs);
}

// Cyclic reduction at the Kokkos team level. Any (thread, vector)
// parameterization is intended to work.
template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void cr (const TeamMember& team,
         TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
         typename std::enable_if<TridiagDiag::rank == 1>::type* = 0,
         typename std::enable_if<DataArray::rank == 1>::type* = 0) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  assert(dl.extent_int(0) == nrow);
  assert(d. extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X. extent_int(0) == nrow);
  const int team_id = get_thread_id_within_team(team);
  const int nteam = get_team_nthr(team);
  int os = 1, stride;
  // Go down reduction.
  while ((stride = (os << 1)) < nrow) {
    const int inc = stride*nteam;
    for (int i = stride*team_id; i < nrow; i += inc) {
      int im = i - os;
      int ip = i + os;
      // GPU does well with ternary ?: op. Use it throughout this
      // impl. It requires the trick noted in a few lines.
      const auto f1 = im >= 0   ? -dl(i)/d(im) : 0;
      const auto f2 = ip < nrow ? -du(i)/d(ip) : 0;
      // Trick to keep im, ip in bounds; the index is modified only
      // when the corresponding f is 0, so the resulting invalid
      // value is multipled by 0.
      im = im >= 0   ? im : i;
      ip = ip < nrow ? ip : i;
      dl(i)  = f1*dl(im);
      du(i)  = f2*du(ip);
      d (i) += f1*du(im) + f2*dl(ip);
      X (i) += f1*X (im) + f2*X (ip);
    }
    os <<= 1;
    // Go down in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
  // Bottom 1 or 2 levels of the reduction. This could be folded into
  // the previous loop, but it's a slight opt to handle these cases
  // separately.
  if (team_id == 0) {
    if (os >= nrow) {
      X(0) /= d(0);
    } else {
      const auto
        det = d(0)*d(os) - du(0)*dl(os),
        x0 = X(0), x1 = X(os);
      X( 0) = (d(os)*x0 - du( 0)*x1)/det;
      X(os) = (d( 0)*x1 - dl(os)*x0)/det;
    }
  }
  team.team_barrier();
  os >>= 1;
  assert(os < nrow);
  // Go up reduction.
  while (os) {
    stride = os << 1;
    const int inc = stride*nteam;
    for (int i = stride*team_id + os; i < nrow; i += inc) {
      const int im = i - os;
      const int ip = i + os;
      assert(im >= 0 || ip < nrow);
      Scalar f = 0;
      f += im >=   0 ? dl(i)*X(im) : 0;
      f += ip < nrow ? du(i)*X(ip) : 0;
      X(i) = (X(i) - f)/d(i);
    }
    os >>= 1;
    // Go up in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void cr (const TeamMember& team,
         TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
         typename std::enable_if<TridiagDiag::rank == 1>::type* = 0,
         typename std::enable_if<DataArray::rank == 2>::type* = 0) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  assert(X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  const int nrhs = X.extent_int(1);
  const int tid = get_thread_id_within_team(team);
  const int nthr = get_team_nthr(team);
  const int team_size = ko::min(nrhs, nthr);
  const int nteam = nthr / team_size;
  const int team_id = tid / team_size;
  const int team_tid = tid % team_size;
  const bool team_lead = tid % team_size == 0;
  int os = 1, stride;
  while ((stride = (os << 1)) < nrow) {
    const int inc = stride*nteam;
    const int nit = (nrow + inc - 1)/inc;
    for (int i = stride*team_id, it = 0; it < nit; i += inc, ++it) {
      int im = i - os;
      int ip = i + os;
      Scalar f1 = 0, f2 = 0;
      const bool run = team_id < nteam && i < nrow;
      assert(team_id != 0 || run);
      if (run) {
        f1 = im >= 0   ? -dl(i)/d(im) : 0;
        f2 = ip < nrow ? -du(i)/d(ip) : 0;
        im = im >= 0   ? im : i;
        ip = ip < nrow ? ip : i;
        for (int j = team_tid; j < nrhs; j += team_size)
          X(i,j) += f1*X(im,j) + f2*X(ip,j);
      }
      // Update A only after all threads are done using current values.
      team.team_barrier();
      if (team_lead && run) {
        dl(i)  = f1*dl(im);
        du(i)  = f2*du(ip);
        d (i) += f1*du(im) + f2*dl(ip);
      }
    }
    os <<= 1;
    team.team_barrier();
  }
  if (team_id == 0) {
    if (os >= nrow) {
      for (int j = team_tid; j < nrhs; j += team_size)
        X(0,j) /= d(0);
    } else {
      for (int j = team_tid; j < nrhs; j += team_size) {
        const auto
          det = d(0)*d(os) - du(0)*dl(os),
          x0 = X(0,j), x1 = X(os,j);
        X( 0,j) = (d(os)*x0 - du( 0)*x1)/det;
        X(os,j) = (d( 0)*x1 - dl(os)*x0)/det;
      }
    }
  }
  team.team_barrier();
  os >>= 1;
  assert(os < nrow);
  while (os) {
    if (team_id < nteam) {
      stride = os << 1;
      const int inc = stride*nteam;
      for (int i = stride*team_id + os; i < nrow; i += inc) {
        const int im = i - os;
        const int ip = i + os;
        assert(im >= 0 || ip < nrow);
        for (int j = team_tid; j < nrhs; j += team_size) {
          Scalar f = 0;
          f += im >=   0 ? dl(i)*X(im,j) : 0;
          f += ip < nrow ? du(i)*X(ip,j) : 0;
          X(i,j) = (X(i,j) - f)/d(i);
        }
      }
    }
    os >>= 1;
    team.team_barrier();
  }
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void cr (const TeamMember& team,
         TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
         typename std::enable_if<TridiagDiag::rank == 2>::type* = 0,
         typename std::enable_if<DataArray::rank == 2>::type* = 0) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(dl.extent_int(1) == nrhs);
  assert(d. extent_int(1) == nrhs);
  assert(du.extent_int(1) == nrhs);
  assert(dl.extent_int(0) == nrow);
  assert(d. extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X. extent_int(0) == nrow);
  const int tid = get_thread_id_within_team(team);
  const int nthr = get_team_nthr(team);
  const int team_size = ko::min(nrhs, nthr);
  const int nteam = nthr / team_size;
  const int team_id = tid / team_size;
  const int team_tid = tid % team_size;
  int os = 1, stride;
  while ((stride = (os << 1)) < nrow) {
    if (team_id < nteam) {
      const int inc = stride*nteam;
      for (int i = stride*team_id; i < nrow; i += inc) {
        int im = i - os;
        int ip = i + os;
        const bool im_ok = im >= 0;
        const bool ip_ok = ip < nrow;
        im = im_ok ? im : i;
        ip = ip_ok ? ip : i;
        for (int j = team_tid; j < nrhs; j += team_size) {
          const auto f1 = im_ok ? -dl(i,j)/d(im,j) : 0;
          const auto f2 = ip_ok ? -du(i,j)/d(ip,j) : 0;
          dl(i,j)  = f1*dl(im,j);
          du(i,j)  = f2*du(ip,j);
          d (i,j) += f1*du(im,j) + f2*dl(ip,j);
          X (i,j) += f1*X (im,j) + f2*X (ip,j);
        }
      }
    }
    os <<= 1;
    team.team_barrier();
  }
  if (team_id == 0) {
    if (os >= nrow) {
      for (int j = team_tid; j < nrhs; j += team_size)
        X(0,j) /= d(0,j);
    } else {
      for (int j = team_tid; j < nrhs; j += team_size) {
        const auto
          det = d(0,j)*d(os,j) - du(0,j)*dl(os,j),
          x0 = X(0,j), x1 = X(os,j);
        X( 0,j) = (d(os,j)*x0 - du( 0,j)*x1)/det;
        X(os,j) = (d( 0,j)*x1 - dl(os,j)*x0)/det;
      }
    }
  }
  team.team_barrier();
  os >>= 1;
  assert(os < nrow);
  while (os) {
    if (team_id < nteam) {
      stride = os << 1;
      const int inc = stride*nteam;
      for (int i = stride*team_id + os; i < nrow; i += inc) {
        const int im = i - os;
        const int ip = i + os;
        assert(im >= 0 || ip < nrow);
        const bool im_ok = im >= 0;
        const bool ip_ok = ip < nrow;
        for (int j = team_tid; j < nrhs; j += team_size) {
          Scalar f = 0;
          f += im_ok ? dl(i,j)*X(im,j) : 0;
          f += ip_ok ? du(i,j)*X(ip,j) : 0;
          X(i,j) = (X(i,j) - f)/d(i,j);
        }
      }
    }
    os >>= 1;
    team.team_barrier();
  }
}

// Pure Cuda version to measure the Kokkos overhead. For this kind of
// computation, every little bit counts, so we need to know that
// overhead.
#ifdef KOKKOS_ENABLE_CUDA
template <typename Scalar>
inline __device__
void cr_a1x1p (const int nrow,
               Scalar* const dl, Scalar* const d, Scalar* const du,
               Scalar* const X) {
  const int team_id = threadIdx.x;
  const int nteam = blockDim.x;
  int os = 1, stride;
  while ((stride = (os << 1)) < nrow) {
    const int inc = stride*nteam;
    for (int i = stride*team_id; i < nrow; i += inc) {
      int im = i - os;
      int ip = i + os;
      const auto f1 = im >= 0   ? -dl[i]/d[im] : 0;
      const auto f2 = ip < nrow ? -du[i]/d[ip] : 0;
      im = im >= 0   ? im : i;
      ip = ip < nrow ? ip : i;
      dl[i]  = f1*dl[im];
      du[i]  = f2*du[ip];
      d [i] += f1*du[im] + f2*dl[ip];
      X [i] += f1*X [im] + f2*X [ip];
    }
    os <<= 1;
    __syncthreads();
  }
  if (team_id == 0) {
    if (os >= nrow) {
      X[0] /= d[0];
    } else {
      const auto
        det = d[0]*d[os] - du[0]*dl[os],
        x0 = X[0], x1 = X[os];
      X[ 0] = (d[os]*x0 - du[ 0]*x1)/det;
      X[os] = (d[ 0]*x1 - dl[os]*x0)/det;
    }
  }
  __syncthreads();
  os >>= 1;
  assert(os < nrow);
  while (os) {
    stride = os << 1;
    for (int i = stride*team_id + os; i < nrow; i += stride*nteam) {
      const int im = i - os;
      const int ip = i + os;
      assert(im >= 0 || ip < nrow);
      Scalar f = 0;
      f += im >= 0   ? dl[i]*X[im] : 0;
      f += ip < nrow ? du[i]*X[ip] : 0;
      X[i] = (X[i] - f)/d[i];
    }
    os >>= 1;
    __syncthreads();
  }
}
#endif // KOKKOS_ENABLE_CUDA

template <typename TridiagDiag, typename XArray, typename YArray>
KOKKOS_INLINE_FUNCTION
int matvec (TridiagDiag dl, TridiagDiag d, TridiagDiag du, XArray X, YArray Y) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X.extent_int(0) == nrow);
  assert(Y.extent_int(0) == nrow);
  assert(Y.extent_int(1) == nrhs);
  if (nrow == 1) {
    for (int j = 0; j < nrhs; ++j)
      Y(0,j) = d(0) * X(0,j);
    return 0;
  }
  for (int j = 0; j < nrhs; ++j)
    Y(0,j) = d(0) * X(0,j) + du(0) * X(1,j);
  for (int i = 1; i < nrow-1; ++i)
    for (int j = 0; j < nrhs; ++j)
      Y(i,j) = (dl(i) * X(i-1,j) +
                d (i) * X(i  ,j) +
                du(i) * X(i+1,j));
  const int i = nrow-1;
  for (int j = 0; j < nrhs; ++j)
    Y(i,j) = dl(i) * X(i-1,j) + d(i) * X(i,j);
  return 0;  
}

KOKKOS_INLINE_FUNCTION
static unsigned int nextpow2 (unsigned int v) {
  // From https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void pcr (const TeamMember& team,
          TridiagDiag dl, TridiagDiag d, TridiagDiag du,
          DataArray X) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  const int nthr = get_team_nthr(team);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X.extent_int(0) == nrow);
  assert(nrow*nrhs <= nthr);
  const int tid = get_thread_id_within_team(team);
  const int team_size = nrhs;
  const int team_id = tid / team_size;
  const int team_tid = tid % team_size;
  const bool team_lead = tid % team_size == 0;
  const int i = team_id;
  const int j = team_tid;
  int n = nextpow2(nrow);
  int os = 1;
  while (n) {
    const int im = i - os;
    const int ip = i + os;
    Scalar a = 0, b, c = 0, f, y;
    if (i < nrow) {
      if (im >= 0 || ip < nrow) {
        b = d(i);
        y = X(i,j);
        if (im >= 0) {
          f  = -dl(i)/d(im);
          a  = f*dl(im);
          b += f*du(im);
          y += f*X (im,j);
        }
        if (ip < nrow) {
          f  = -du(i)/d(ip);
          b += f*dl(ip);
          c  = f*du(ip);
          y += f*X (ip,j);
        }
      } else {      
        b = 1;
        y = X(i,j)/d(i);
      }
    }
    team.team_barrier();
    if (i < nrow) {
      X(i,j) = y;
      if (team_lead) {
        dl(i) = a;
        d (i) = b;
        du(i) = c;
      }
    }
    team.team_barrier();
    os <<= 1;
    n >>= 1;
  }
}

#ifdef TRIDIAG_HAVE_CUSPARSE
template <typename T>
cusparseStatus_t cusparseTgtsv2StridedBatch_bufferSizeExt(
  cusparseHandle_t, int, const T*, const T*, const T*, const T*, int, int, size_t*);
template <typename T>
cusparseStatus_t cusparseTgtsv2StridedBatch(
  cusparseHandle_t, int, const T*, const T*, const T*, T*, int, int, T*);
template <> cusparseStatus_t cusparseTgtsv2StridedBatch_bufferSizeExt<float> (
  cusparseHandle_t h, int m, const float* dl, const float* d, const float* du, const float* x,
  int bc, int bs, size_t* sz)
{ return cusparseSgtsv2StridedBatch_bufferSizeExt(h, m, dl, d, du, x, bc, bs, sz); }
template <> cusparseStatus_t cusparseTgtsv2StridedBatch<float> (
  cusparseHandle_t h, int m, const float* dl, const float* d, const float* du, float* x,
  int bc, int bs, float* b)
{ return cusparseSgtsv2StridedBatch(h, m, dl, d, du, x, bc, bs, b); }
template <> cusparseStatus_t cusparseTgtsv2StridedBatch_bufferSizeExt<double> (
  cusparseHandle_t h, int m, const double* dl, const double* d, const double* du, const double* x,
  int bc, int bs, size_t* sz)
{ return cusparseDgtsv2StridedBatch_bufferSizeExt(h, m, dl, d, du, x, bc, bs, sz); }
template <> cusparseStatus_t cusparseTgtsv2StridedBatch<double> (
  cusparseHandle_t h, int m, const double* dl, const double* d, const double* du, double* x,
  int bc, int bs, double* b)
{ return cusparseDgtsv2StridedBatch(h, m, dl, d, du, x, bc, bs, b); }

template <typename Scalar>
class CusparseSolver {
  int nrow_, nprob_;
  bool ok_;
  cusparseHandle_t handle_;
  Kokkos::View<Scalar*> buffer_;

  void init (const int& nrow, const int& nprob) {
    nrow_ = nrow;
    nprob_ = nprob;
    ok_ = true;
    if (nrow_ < 3) ok_ = false;
    auto status = cusparseCreate(&handle_);
    if (status != CUSPARSE_STATUS_SUCCESS) {
      ok_ = false;
      return;
    }
    size_t sz;
    status = cusparseTgtsv2StridedBatch_bufferSizeExt<Scalar>(
      handle_, nrow_, nullptr, nullptr, nullptr, nullptr, nprob_, nrow_, &sz);
    if (status != CUSPARSE_STATUS_SUCCESS) {
      ok_ = false;
      return;
    }
    buffer_ = Kokkos::View<Scalar*>("buffer", sz/sizeof(Scalar));
  }

public:
  CusparseSolver (const int& nrow, const int& nprob) {
    init(nrow, nprob);
    assert(ok_);
  }

  ~CusparseSolver () {
    cusparseDestroy(handle_);
  }

  template <typename Acs, typename Dcs>
  bool solve (const Acs& A, Dcs& X) {
    if ( ! ok_) return false;
    auto status = cusparseTgtsv2StridedBatch(
      handle_, nrow_, &A.impl_map().reference(0,0,0), &A.impl_map().reference(1,0,0),
      &A.impl_map().reference(2,0,0), X.data(), nprob_, nrow_, buffer_.data());
    return status == CUSPARSE_STATUS_SUCCESS;
  }

  template <typename Astd, typename Acs>
  void to_matrix (const Astd& a, const Acs& b, const int nrhs) const {
    const int nprob = nprob_;
    const int nrow = nrow_;
    const int nrp = nrow*nprob;
    assert(a.extent_int(0) == nprob/nrhs);
    assert(a.extent_int(2) == nrow);
    assert(a.extent_int(1) == 3);
    assert(b.extent_int(0) == 3);
    assert(b.extent_int(1) == nprob);
    assert(b.extent_int(2) == nrow);
    auto f = KOKKOS_LAMBDA (const int i) {
      const int c = i/nrp;
      const int prob = (i % nrp) / nrow;
      const int row = i % nrow;
      b(c,prob,row) = a(prob/nrhs,c,row);
    };
    Kokkos::parallel_for(3*nrow*nprob, f);
  }

  template <typename Dstd, typename Dcs>
  void to_data (const Dstd& a, const Dcs& b) const {
    const int nrhs = a.extent_int(2);
    const int nprob = nprob_;
    const int nrow = nrow_;
    assert(a.extent_int(0) == nprob/nrhs);
    assert(a.extent_int(1) == nrow);
    assert(a.extent_int(2) == nrhs);
    assert(b.extent_int(0) == nprob);
    assert(b.extent_int(1) == nrow);
    auto f = KOKKOS_LAMBDA (const int i) {
      const int prob = i / nrow;
      const int row = i % nrow;
      const int rhs = prob % nrhs;
      b(prob,row) = a(prob/nrhs,row,rhs);
    };
    Kokkos::parallel_for(nrow*nprob, f);
  }

  template <typename Dcs, typename Dstd>
  void from_data (const Dcs& a, const Dstd& b) const {
    const int nrhs = b.extent_int(2);
    const int nprob = nprob_;
    const int nrow = nrow_;
    assert(b.extent_int(0) == nprob/nrhs);
    assert(b.extent_int(1) == nrow);
    assert(b.extent_int(2) == nrhs);
    assert(a.extent_int(0) == nprob);
    assert(a.extent_int(1) == nrow);
    auto f = KOKKOS_LAMBDA (const int i) {
      const int prob = i / nrow;
      const int row = i % nrow;
      const int rhs = prob % nrhs;
      b(prob/nrhs,row,rhs) = a(prob,row);
    };
    Kokkos::parallel_for(nrow*nprob, f);
  }
};
#endif // TRIDIAG_HAVE_CUSPARSE

#if defined TRIDIAG_HAVE_LAPACK || defined TRIDIAG_HAVE_MKL
# ifndef TRIDIAG_HAVE_MKL
extern "C" {
void sgttrf_(int*, float*, float*, float*, float*, int*, int*);
void sgttrs_(char*, int*, int*, float*, float*, float*, float*, int*, float*, int*, int*);
void dgttrf_(int*, double*, double*, double*, double*, int*, int*);
void dgttrs_(char*, int*, int*, double*, double*, double*, double*, int*, double*, int*, int*);
}
# endif // TRIDIAG_HAVE_MKL

template <typename T>
int gttrf(int n, T* dl, T* d, T* du, T* du2, int* ipiv);
template <typename T>
int gttrs(char dir, int n, int nrhs, T* dl, T* d, T* du, T* du2, int* ipiv,
          T* x, int ldx);

template <>
inline int gttrf<float> (int n, float* dl, float* d, float* du, float* du2, int* ipiv) {
  int info;
  sgttrf_(&n, dl, d, du, du2, ipiv, &info);
  return info;
}
template <>
inline int gttrs<float> (char dir, int n, int nrhs, float* dl, float* d, float* du,
                         float* du2, int* ipiv, float* x, int ldx) {
  int info;
  sgttrs_(&dir, &n, &nrhs, dl, d, du, du2, ipiv, x, &ldx, &info);
  return info;
}

template <>
inline int gttrf<double> (int n, double* dl, double* d, double* du, double* du2, int* ipiv) {
  int info;
  dgttrf_(&n, dl, d, du, du2, ipiv, &info);
  return info;
}
template <>
inline int gttrs<double> (char dir, int n, int nrhs, double* dl, double* d, double* du,
                          double* du2, int* ipiv, double* x, int ldx) {
  int info;
  dgttrs_(&dir, &n, &nrhs, dl, d, du, du2, ipiv, x, &ldx, &info);
  return info;
}
#endif // TRIDIAG_HAVE_LAPACK || defined TRIDIAG_HAVE_MKL

#if defined TRIDIAG_HAVE_MKL
# ifndef TRIDIAG_HAVE_MKL
extern "C" {
void sdttrfb_(int*, float*, float*, float*, int*);
void sdttrsb_(char*, int*, int*, float*, float*, float*, float*, int*, int*);
void ddttrfb_(int*, double*, double*, double*, int*);
void ddttrsb_(char*, int*, int*, double*, double*, double*, double*, int*, int*);
}
#endif // TRIDIAG_HAVE_MKL

template <typename T>
int dttrfb(int n, T* dl, T* d, T* du);
template <typename T>
int dttrsb(char dir, int n, int nrhs, T* dl, T* d, T* du, T* x, int ldx);

template <>
inline int dttrfb<float> (int n, float* dl, float* d, float* du) {
  int info;
  sdttrfb_(&n, dl, d, du, &info);
  return info;
}
template <>
inline int dttrsb<float> (char dir, int n, int nrhs, float* dl, float* d, float* du,
                          float* x, int ldx) {
  int info;
  sdttrsb_(&dir, &n, &nrhs, dl, d, du, x, &ldx, &info);
  return info;
}

template <>
inline int dttrfb<double> (int n, double* dl, double* d, double* du) {
  int info;
  ddttrfb_(&n, dl, d, du, &info);
  return info;
}
template <>
inline int dttrsb<double> (char dir, int n, int nrhs, double* dl, double* d, double* du,
                           double* x, int ldx) {
  int info;
  ddttrsb_(&dir, &n, &nrhs, dl, d, du, x, &ldx, &info);
  return info;
}
#endif // TRIDIAG_HAVE_MKL

template <typename TridiagDiag>
KOKKOS_INLINE_FUNCTION
void fill_tridiag_matrix (TridiagDiag dl, TridiagDiag d, TridiagDiag du,
                          const int& seed) {
  const int nrow = d.extent_int(0);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  for (int i = 0; i < nrow; ++i) {
    const int k = seed + i;
    dl(i) = (k % 5 == 0 ? -1 : 1) * 1.3 * (0.1 + ((k*k) % 11));
    du(i) = (k % 7 == 0 ? -1 : 1) * 1.7 * (0.2 + ((k*k) % 13));
    d (i) = ((k % 3 == 0 ? -1 : 1) *
             (0.7 + std::abs(dl(i)) + std::abs(du(i)) + (k % 17)));
  }
}

template <typename DataArray>
KOKKOS_INLINE_FUNCTION
void fill_data_matrix (DataArray X, const int& seed) {
  const int nrow = X.extent_int(0);
  const int nrhs = X.extent_int(1);
  for (int i = 0; i < nrow; ++i)
    for (int j = 0; j < nrhs; ++j)
      X(i,j) = (((7*i + 11*j + 3*i*j) % 3 == 0 ? -1 : 1) *
                1.7 * ((17*(i - 19) + 13*(j - 11) + 5*(i - 5)*(j - 7) + seed) % 47));
}

#ifndef TRIDIAG_DOUBLE_PRECISION
# define TRIDIAG_DOUBLE_PRECISION 1
#endif
#if TRIDIAG_DOUBLE_PRECISION
typedef double Real;
#else
typedef float Real;
#endif

template <typename Array>
Real reldif (const Array& a, const Array& b) {
  assert(a.extent_int(0) == b.extent_int(0));
  assert(a.extent_int(1) == b.extent_int(1));
  assert(a.rank == 2);
  assert(b.rank == 2);
  Real num = 0, den = 0;
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j) {
      num = std::max(num, std::abs(a(i,j) - b(i,j)));
      den = std::max(den, std::abs(a(i,j)));
    }
  return num/den;
}

struct Solver {
  enum Enum { thomas, thomas_pack_a1xm, thomas_pack_amxm, thomas_amxm,
              gttr, dttr,
              cr_a1xm, cr_amxm, cr_a1x1, cr_a1x1p, pcr, cusparse,
              error };
  static std::string convert (Enum e) {
    switch (e) {
    case thomas: return "thomas";
    case thomas_pack_a1xm: return "thomas_pack_a1xm";
    case thomas_pack_amxm: return "thomas_pack_amxm";
    case thomas_amxm: return "thomas_amxm";
    case gttr: return "gttr";
    case dttr: return "dttr";
    case cr_a1xm: return "cr_a1xm";
    case cr_amxm: return "cr_amxm";
    case cr_a1x1: return "cr_a1x1";
    case cr_a1x1p: return "cr_a1x1p";
    case pcr: return "pcr";
    case cusparse: return "cusparse";
    default:
      assert(0);
      return "";
    }
  }
  static Enum convert (const std::string& s) {
    if (s == "thomas") return thomas;
    if (s == "thomas_pack_a1xm") return thomas_pack_a1xm;
    if (s == "thomas_pack_amxm") return thomas_pack_amxm;
    if (s == "thomas_amxm") return thomas_amxm;
    if (s == "gttr") return gttr;
    if (s == "dttr") return dttr;
    if (s == "cr_a1xm") return cr_a1xm;
    if (s == "cr_amxm") return cr_amxm;
    if (s == "cr_a1x1") return cr_a1x1;
    if (s == "cr_a1x1p") return cr_a1x1p;
    if (s == "pcr") return pcr;
    if (s == "cusparse") return cusparse;
    return error;
  }
  static Enum all[];
};

Solver::Enum Solver::all[] = {
  thomas, cr_a1xm, cr_amxm, cr_a1x1, pcr
#ifdef TRIDIAG_HAVE_CUSPARSE
  //, cusparse
#endif
};

// LayoutStride has too much overhead to use it. Instead we use the
// (dl,d,du) triple rather than a matrix A to get the layout
// flexibility we need.
//#define TRIDIAG_STRIDE
#ifdef TRIDIAG_STRIDE
using BulkLayout = Kokkos::LayoutStride;
using TeamLayout = Kokkos::LayoutLeft;
#else
using BulkLayout = Kokkos::LayoutRight;
using TeamLayout = Kokkos::LayoutRight;
#endif

template <typename ScalarType>
using TridiagArray = Kokkos::View<ScalarType*[3], TeamLayout>;
template <typename ScalarType>
using DataArray = Kokkos::View<ScalarType**, TeamLayout>;

static void test1 () {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using MT = typename TeamPolicy::member_type;
  const auto policy = get_default_team_policy<>(1, 128, 1);
  const int nthr = policy.team_size();
  for (const auto solver : Solver::all) {
    for (int nrow = 1; nrow <= 129; ++nrow) {
      const int nrhs_max = solver != Solver::pcr ? 60 : 8;
      const int nrhs_inc = solver != Solver::pcr ? 11 : 1;
      for (int nrhs = 1; nrhs <= nrhs_max; nrhs += nrhs_inc) {
        if (solver == Solver::pcr && nrhs*nrow >= nthr) continue;
        if (solver == Solver::cr_amxm && nrhs > 1) continue;
        if (solver == Solver::cr_a1x1 && nrhs > 1) continue;
        TridiagArray<Real> A("A", nrow), Acopy("A", A.extent(0));
        DataArray<Real> B("B", A.extent(0), nrhs), X("X", B.extent(0), B.extent(1)),
          Y("Y", X.extent(0), X.extent(1));
        auto Am = Kokkos::create_mirror_view(A);
        auto Bm = Kokkos::create_mirror_view(B);
        {
          const auto dl = Kokkos::subview(Am, Kokkos::ALL(), 0);
          const auto d  = Kokkos::subview(Am, Kokkos::ALL(), 1);
          const auto du = Kokkos::subview(Am, Kokkos::ALL(), 2);
          fill_tridiag_matrix(dl, d, du, nrhs);
        }
        fill_data_matrix(Bm, nrhs);
        Kokkos::deep_copy(A, Am);
        Kokkos::deep_copy(B, Bm);
        Kokkos::deep_copy(Acopy, A);
        Kokkos::deep_copy(X, B);
        const auto policy = get_default_team_policy<>(1, nrow, nrhs);
        switch (solver) {
        case Solver::thomas: {
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            const auto dl = Kokkos::subview(A, Kokkos::ALL(), 0);
            const auto d  = Kokkos::subview(A, Kokkos::ALL(), 1);
            const auto du = Kokkos::subview(A, Kokkos::ALL(), 2);
            thomas(team, dl, d, du, X);
          };
          Kokkos::parallel_for(policy, f);
        } break;
        case Solver::cr_a1xm: {
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            const auto dl = Kokkos::subview(A, Kokkos::ALL(), 0);
            const auto d  = Kokkos::subview(A, Kokkos::ALL(), 1);
            const auto du = Kokkos::subview(A, Kokkos::ALL(), 2);
            cr(team, dl, d, du, X);
          };
          Kokkos::parallel_for(policy, f);
        } break;
        case Solver::cr_amxm: {
          Kokkos::View<Real***, TeamLayout> Am(A.data(), nrow, 3, nrhs);
          const auto dl = Kokkos::subview(Am, Kokkos::ALL(), 0, Kokkos::ALL());
          const auto d  = Kokkos::subview(Am, Kokkos::ALL(), 1, Kokkos::ALL());
          const auto du = Kokkos::subview(Am, Kokkos::ALL(), 2, Kokkos::ALL());
          const auto f = KOKKOS_LAMBDA (const MT& team) { cr(team, dl, d, du, X); };
          Kokkos::parallel_for(policy, f);
        } break;
        case Solver::cr_a1x1: {
          Kokkos::View<Real**, TeamLayout> Am(A.data(), nrow, 3);
          const auto dl = Kokkos::subview(Am, Kokkos::ALL(), 0);
          const auto d  = Kokkos::subview(Am, Kokkos::ALL(), 1);
          const auto du = Kokkos::subview(Am, Kokkos::ALL(), 2);
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            cr(team, dl, d, du, Kokkos::subview(X, Kokkos::ALL(), 0));
          };
          Kokkos::parallel_for(policy, f);
        } break;
        case Solver::pcr: {
          Kokkos::View<Real**, TeamLayout> Am(A.data(), nrow, 3);
          const auto dl = Kokkos::subview(Am, Kokkos::ALL(), 0);
          const auto d  = Kokkos::subview(Am, Kokkos::ALL(), 1);
          const auto du = Kokkos::subview(Am, Kokkos::ALL(), 2);
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            pcr(team, dl, d, du, X);
          };
          Kokkos::parallel_for(policy, f);
        } break;
        default:
          std::cout << "test1 does not support " << Solver::convert(solver)
                    << "\n";
        }
        Real re; {
          auto Acopym = Kokkos::create_mirror_view(Acopy);
          const auto dl = Kokkos::subview(Acopym, Kokkos::ALL(), 0);
          const auto d  = Kokkos::subview(Acopym, Kokkos::ALL(), 1);
          const auto du = Kokkos::subview(Acopym, Kokkos::ALL(), 2);
          auto Xm = Kokkos::create_mirror_view(X);
          auto Ym = Kokkos::create_mirror_view(Y);
          Kokkos::deep_copy(Acopym, Acopy);
          Kokkos::deep_copy(Xm, X);
          matvec(dl, d, du, Xm, Ym);
          re = reldif(Bm, Ym);
        }
        if (re > 50*std::numeric_limits<Real>::epsilon())
          std::cout << "test1: solver " << Solver::convert(solver)
                    << " nrow " << nrow << " nrhs " << nrhs
                    << " re " << re << "\n";
      }
    }
  }
}

struct Input {
  Solver::Enum method;
  int nprob, nrow, nrhs, nwarp;
  bool use_scratch, gpu;

  Input ()
    : method(Solver::thomas), nprob(2048), nrow(128), nrhs(43), nwarp(-1),
      use_scratch(false), gpu(true)
  {}

  bool parse (int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
      if (util::eq(argv[i], "-m", "--method")) {
        util::expect_another_arg(i, argc);
        method = Solver::convert(argv[++i]);
        if (method == Solver::error) {
          std::cout << "Not a solver: " << argv[i] << "\n";
          return false;
        }
      } else if (util::eq(argv[i], "-np", "--nprob")) {
        util::expect_another_arg(i, argc);
        nprob = std::atoi(argv[++i]);
      } else if (util::eq(argv[i], "-nr", "--nrow")) {
        util::expect_another_arg(i, argc);
        nrow = std::atoi(argv[++i]);
      } else if (util::eq(argv[i], "-nc", "--nrhs")) {
        util::expect_another_arg(i, argc);
        nrhs = std::atoi(argv[++i]);
      } else if (util::eq(argv[i], "-nw", "--nwarp")) {
        util::expect_another_arg(i, argc);
        nwarp = std::atoi(argv[++i]);
      } else if (util::eq(argv[i], "-ng", "--notgpu")) {
        gpu = false;
      } else {
        std::cout << "Unexpected arg: " << argv[i] << "\n";
        return false;
      }
    }
    return true;
  }
};

std::string string (const Input& in) {
  std::stringstream ss;
  ss << "run: solver " << Solver::convert(in.method) << " nprob " << in.nprob
     << " nrow " << in.nrow << " nrhs " << in.nrhs << "\n";
  return ss.str();
}

#ifdef KOKKOS_ENABLE_CUDA
template <typename Real>
__global__
void call_cr_a1x1p (const int nrow, Real* const Adata, Real* const Xdata) {
  const int ip = blockIdx.x;
  Real *dl = Adata + ip*3*nrow, *d = dl+nrow, *du = d+nrow, *x = Xdata + ip*nrow;
  cr_a1x1p(nrow, dl, d, du, x);
}

template <typename AT, typename XT>
double dispatch_cr_a1x1p (const Input& in, AT& A, XT& X) {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  TeamPolicy policy(in.nprob, 32*in.nwarp, 1);
  Kokkos::fence();
  double t0 = util::gettime();
  const int nrow = in.nrow;
  Real* const Adata = A.data();
  Real* const Xdata = X.data();
  cudaFuncSetCacheConfig(call_cr_a1x1p<Real>, cudaFuncCachePreferL1);
  call_cr_a1x1p<<<in.nprob,32*in.nwarp>>>(nrow, Adata, Xdata);
  Kokkos::fence();
  double t1 = util::gettime();
  return t1 - t0;
}
#endif // KOKKOS_ENABLE_CUDA

template <typename ScalarType>
using TridiagArrays = Kokkos::View<ScalarType***, BulkLayout>;
template <typename ScalarType>
using DataArrays = Kokkos::View<ScalarType***, BulkLayout>;
template <typename ScalarType>
using TridiagArraysUm = Kokkos::View<ScalarType***, BulkLayout, Kokkos::MemoryUnmanaged>;
template <typename ScalarType>
using DataArraysUm = Kokkos::View<ScalarType***, BulkLayout, Kokkos::MemoryUnmanaged>;

static void run_gpu (const Input& in) {
  using Kokkos::subview;
  using Kokkos::ALL;
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  std::cout << string(in);
#ifdef TRIDIAG_STRIDE
  const int order[] = {2,0,1};
  const int Adim[] = {in.nprob, in.nrow, 3};
  const int Xdim[] = {in.nprob, in.nrow, in.nrhs};
  const auto Astride = Kokkos::LayoutStride::order_dimensions(3, order, Adim);
  const auto Xstride = Kokkos::LayoutStride::order_dimensions(3, order, Xdim);
  TridiagArrays<Real> A("A", Astride), Acopy("Acopy", Astride);
  DataArrays<Real> B("B", Xstride), X("X", Xstride), Y("Y", Xstride);
#else
  TridiagArrays<Real> A("A", in.nprob, 3, in.nrow),
    Acopy("Acopy", A.extent(0), A.extent(1), A.extent(2));
  DataArrays<Real> B("B", in.nprob, in.nrow, in.nrhs), X("X", in.nprob, in.nrow, in.nrhs),
    Y("Y", in.nprob, in.nrow, in.nrhs);
#endif
  auto Am = Kokkos::create_mirror_view(A);
  auto Bm = Kokkos::create_mirror_view(B);
  for (int i = 0; i < in.nprob; ++i) {
    const auto dl = subview(Am, i, 0, ALL());
    const auto d  = subview(Am, i, 1, ALL());
    const auto du = subview(Am, i, 2, ALL());
    fill_tridiag_matrix(dl, d, du, i);
    fill_data_matrix(subview(Bm, i, ALL(), ALL()), i);
  }
  Kokkos::deep_copy(A, Am);
  Kokkos::deep_copy(B, Bm);
  Kokkos::deep_copy(Acopy, A);
  Kokkos::deep_copy(X, B);
  auto policy = in.nwarp < 0 ? get_default_team_policy<>(in.nprob, in.nrow, in.nrhs) :
    TeamPolicy(in.nprob, 32*in.nwarp, 1);
  if (in.nwarp > 0) {
    assert(policy.team_size() == 32*in.nwarp);
  } else {
    std::cout << "run_gpu: league " << policy.league_size() << " team "
              << policy.team_size() << "\n";
  }
  if (in.method == Solver::pcr && policy.team_size() < in.nrow*in.nrhs) {
    std::stringstream ss;
    ss << "PCR requires nthr >= nrow*nrhs but nthr " << policy.team_size()
       << " nrow " << in.nrow << " nrhs " << in.nrhs << "\n";
    throw std::runtime_error(ss.str());
  }
  if (in.method == Solver::cr_a1x1 && in.nrhs != 1)
    throw std::runtime_error("In this performance test, cr_a1x1 require nrhs 1.");
  Kokkos::fence();
  double t0, t1;
  switch (in.method) {
  case Solver::thomas: {
    t0 = util::gettime();
    const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
      const int ip = team.league_rank();
      const auto dl = subview(A, ip, 0, ALL());
      const auto d  = subview(A, ip, 1, ALL());
      const auto du = subview(A, ip, 2, ALL());
      thomas(team, dl, d, du, subview(X, ip, ALL(), ALL()));
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    t1 = util::gettime();
  } break;
  case Solver::cr_a1xm: {
    t0 = util::gettime();
    const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
      const int ip = team.league_rank();
      const auto dl = subview(A, ip, 0, ALL());
      const auto d  = subview(A, ip, 1, ALL());
      const auto du = subview(A, ip, 2, ALL());
      cr(team, dl, d, du, subview(X, ip, ALL(), ALL()));
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    t1 = util::gettime();
  } break;
  case Solver::cr_amxm: {
    if (in.nrhs == 1) {
      t0 = util::gettime();
      const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
        const int ip = team.league_rank();
        Kokkos::View<Real***, TeamLayout, Kokkos::MemoryUnmanaged>
          Am(&A(ip,0,0), 3, in.nrow, 1);
        const auto dl = subview(Am, 0, ALL(), ALL());
        const auto d  = subview(Am, 1, ALL(), ALL());
        const auto du = subview(Am, 2, ALL(), ALL());
        cr(team, dl, d, du, subview(X, ip, ALL(), ALL()));
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      t1 = util::gettime();
    } else {
      Kokkos::View<Real****, Kokkos::LayoutRight> Arep("Arep", in.nprob, 3, in.nrow, in.nrhs);
      const int nprob = in.nprob, nrow = in.nrow, nrhs = in.nrhs;
      const int ncrr = 3*nrow*nrhs, nrr = nrow*nrhs;
      auto convert = KOKKOS_LAMBDA (const int i) {
        const int prob = i/ncrr;
        const int c = (i % ncrr) / nrr;
        const int row = (i % nrr) / nrhs;
        const int rhs = i % nrhs;
        Arep(prob,c,row,rhs) = A(prob,c,row);
      };
      Kokkos::parallel_for(nprob*3*nrow*nrhs, convert);
      Kokkos::fence();
      t0 = util::gettime();
      const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
        const int ip = team.league_rank();
        const auto dl = subview(Arep, ip, 0, ALL(), ALL());
        const auto d  = subview(Arep, ip, 1, ALL(), ALL());
        const auto du = subview(Arep, ip, 2, ALL(), ALL());
        cr(team, dl, d, du, subview(X, ip, ALL(), ALL()));
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      t1 = util::gettime();      
    }
  } break;
  case Solver::cr_a1x1: {
    t0 = util::gettime();
    const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
      const int ip = team.league_rank();
      const auto dl = subview(A, ip, 0, ALL());
      const auto d  = subview(A, ip, 1, ALL());
      const auto du = subview(A, ip, 2, ALL());
      const auto x = subview(X, ip, ALL(), 0);
      cr(team, dl, d, du, x);
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    t1 = util::gettime();
  } break;
#ifdef KOKKOS_ENABLE_CUDA
  case Solver::cr_a1x1p: {
    t0 = 0;
    t1 = dispatch_cr_a1x1p(in, A, X);
  } break;
#endif
  case Solver::pcr: {
    t0 = util::gettime();
    const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
      const int ip = team.league_rank();
      const auto dl = subview(A, ip, 0, ALL());
      const auto d  = subview(A, ip, 1, ALL());
      const auto du = subview(A, ip, 2, ALL());
      pcr(team, dl, d, du, subview(X, ip, ALL(), ALL()));
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    t1 = util::gettime();
  } break;
#ifdef TRIDIAG_HAVE_CUSPARSE
  case Solver::cusparse: {
    Kokkos::View<Real***, Kokkos::LayoutRight> Acs("Acs", 3, in.nrhs*in.nprob, in.nrow);
    Kokkos::View<Real**, Kokkos::LayoutRight> Xcs("Xcs", in.nrhs*in.nprob, in.nrow);
    CusparseSolver<Real> cs(in.nrow, in.nrhs*in.nprob);
    cs.to_matrix(A, Acs, in.nrhs);
    cs.to_data(X, Xcs);
    Kokkos::fence();
    t0 = util::gettime();
    cs.solve(Acs, Xcs);
    Kokkos::fence();
    t1 = util::gettime();
    cs.from_data(Xcs, X);
  } break;
#endif
  default:
    std::cout << "run_gpu does not support " << Solver::convert(in.method)
              << "\n";
  }
  const auto et = t1 - t0;
  printf("run: et %1.3e et/datum %1.3e\n", et, et/(in.nprob*in.nrow*in.nrhs));
  Real re; {
    auto Acopym = Kokkos::create_mirror_view(Acopy);
    auto Xm = Kokkos::create_mirror_view(X);
    auto Ym = Kokkos::create_mirror_view(Y);
    Kokkos::deep_copy(Acopym, Acopy);
    Kokkos::deep_copy(Xm, X);
    const auto ip = std::max(0, in.nprob-1);
    const auto dl = subview(Acopym, ip, 0, ALL());
    const auto d  = subview(Acopym, ip, 1, ALL());
    const auto du = subview(Acopym, ip, 2, ALL());
    matvec(dl, d, du,
           subview(Xm, in.nprob-1, ALL(), ALL()),
           subview(Ym, in.nprob-1, ALL(), ALL()));
    re = reldif(
      subview(Bm, in.nprob-1, ALL(), ALL()),
      subview(Ym, in.nprob-1, ALL(), ALL()));
  }
  if (re > 50*std::numeric_limits<Real>::epsilon())
    std::cout << "run: " << " re " << re << "\n";
}

template <typename ST, typename DT>
void transpose_data (const ST& a, DT b) {
  const int np = a.extent_int(0);
  const int m = a.extent_int(1), n = a.extent_int(2), mn = m*n;
  assert(b.extent_int(0) == np);
  assert(b.extent_int(1) == n);
  assert(b.extent_int(2) == m);
  auto f = KOKKOS_LAMBDA (const int i) {
    const int prob = i / mn;
    const int row = (i % mn) / n;
    const int col = i % n;
    b(prob,col,row) = a(prob,row,col);
  };
  Kokkos::parallel_for(np*m*n, f);
}

template <typename ST, typename DT>
void pack_data (const ST& a, DT b) {
  using Pack = typename DT::non_const_value_type;
  static_assert(Pack::packtag, "DT value type must be Pack");
  static_assert(std::is_same<typename ST::non_const_value_type,
                             typename Pack::scalar>::value,
                "ST value and DT::Pack::scalar types must be the same");
  const int np = a.extent_int(0);
  const int m = a.extent_int(1), n = a.extent_int(2), mn = m*n;
  assert(b.extent_int(0) == np);
  assert(b.extent_int(1) == m);
  assert(b.extent_int(2) == (n + Pack::n - 1)/Pack::n);
  auto f = KOKKOS_LAMBDA (const int i) {
    const int prob = i / mn;
    const int row = (i % mn) / n;
    const int col = i % n;
    b(prob, row, col / Pack::n)[col % Pack::n] = a(prob,row,col);
  };
  Kokkos::parallel_for(np*m*n, f);
}

template <typename ST, typename DT>
void unpack_data (const ST& a, DT b) {
  using Pack = typename ST::non_const_value_type;
  static_assert(Pack::packtag, "ST value type must be Pack");
  static_assert(std::is_same<typename DT::non_const_value_type,
                             typename Pack::scalar>::value,
                "DT value and ST::Pack::scalar types must be the same");
  const int np = b.extent_int(0);
  const int m = b.extent_int(1), n = b.extent_int(2), mn = m*n;
  assert(a.extent_int(0) == np);
  assert(a.extent_int(1) == m);
  assert(a.extent_int(2) == (n + Pack::n - 1)/Pack::n);
  auto f = KOKKOS_LAMBDA (const int i) {
    const int prob = i / mn;
    const int row = (i % mn) / n;
    const int col = i % n;
    b(prob,row,col) = a(prob, row, col / Pack::n)[col % Pack::n];
  };
  Kokkos::parallel_for(np*m*n, f);
}

template <typename ST, typename DT>
void pack_matrix (const ST& a, DT b, const int nrhs) {
  using Pack = typename DT::non_const_value_type;
  static_assert(Pack::packtag, "DT value type must be Pack");
  static_assert(std::is_same<typename ST::non_const_value_type,
                             typename Pack::scalar>::value,
                "ST value and DT::Pack::scalar types must be the same");
  const int np = a.extent_int(0);
  const int nrow = a.extent_int(2);
  assert(b.extent_int(0) == np);
  assert(b.extent_int(1) == 3);
  assert(b.extent_int(2) == nrow);
  assert(b.extent_int(3) == (nrhs + Pack::n - 1)/Pack::n);
  auto f = KOKKOS_LAMBDA (const int i) {
    const int prob = i / (3*nrow*nrhs);
    const int c = (i % (3*nrow*nrhs)) / (nrow*nrhs);
    const int row = (i % (nrow*nrhs)) / nrhs;
    const int rhs = i % nrhs;
    b(prob, c, row, rhs / Pack::n)[rhs % Pack::n] = a(prob,c,row);
  };
  Kokkos::parallel_for(np*3*nrow*nrhs, f);
}

template <typename ST, typename DT>
void pack_scalar_matrix (const ST& a, DT b, const int nrhs) {
  const int np = a.extent_int(0);
  const int nrow = a.extent_int(2);
  assert(b.extent_int(0) == np);
  assert(b.extent_int(1) == 3);
  assert(b.extent_int(2) == nrow);
  assert(b.extent_int(3) == nrhs);
  auto f = KOKKOS_LAMBDA (const int i) {
    const int prob = i / (3*nrow*nrhs);
    const int c = (i % (3*nrow*nrhs)) / (nrow*nrhs);
    const int row = (i % (nrow*nrhs)) / nrhs;
    const int rhs = i % nrhs;
    b(prob,c,row,rhs) = a(prob,c,row);
  };
  Kokkos::parallel_for(np*3*nrow*nrhs, f);
}

static void run_cpu (const Input& in) {
  using Pack = pack::Pack<Real,8>;
  using Kokkos::subview;
  using Kokkos::ALL;
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  std::cout << string(in);
  TridiagArrays<Real> A("A", in.nprob, 3, in.nrow),
    Acopy("Acopy", A.extent(0), A.extent(1), A.extent(2));
  DataArrays<Real> B("B", in.nprob, in.nrow, in.nrhs), Y("Y", in.nprob, in.nrow, in.nrhs),
    X("X", in.nprob, in.nrow, in.nrhs);
  Kokkos::View<Real***, BulkLayout> Xt("Xt", in.nprob, in.nrhs, in.nrow);
# ifdef KOKKOS_ENABLE_OPENMP
# pragma omp parallel for
# endif
  for (int i = 0; i < in.nprob; ++i) {
    const auto dl = subview(A, i, 0, ALL());
    const auto d  = subview(A, i, 1, ALL());
    const auto du = subview(A, i, 2, ALL());
    fill_tridiag_matrix(dl, d, du, i);
    fill_data_matrix(subview(X, i, ALL(), ALL()), i);
  }
  Kokkos::deep_copy(Acopy, A);
  Kokkos::deep_copy(B, X);
  Kokkos::fence();
  const auto policy = get_default_team_policy<>(in.nprob, in.nrow, in.nrhs);
#ifdef TRIDIAG_HAVE_MKL
  mkl_set_num_threads(1);
#endif
  double t0, t1;
  switch (in.method) {
  case Solver::thomas: {
    for (int trial = 0; trial < 2; ++trial) {
      Kokkos::deep_copy(A, Acopy);
      // Kokkos::deep_copy was messing up the time for some reason. Do
      // it manually.
      const auto dc = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
          const int i = team.league_rank();
          for (int r = 0; r < in.nrow; ++r)
            for (int c = 0; c < in.nrhs; ++c)
              X(i,r,c) = B(i,r,c);
        });
      };
      Kokkos::parallel_for(policy, dc);
      Kokkos::fence();
      t0 = util::gettime();
      if (in.nrhs == 1) {
        const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
          Kokkos::single(Kokkos::PerTeam(team), [&] () {
            const int ip = team.league_rank();
            const auto dl = subview(A, ip, 0, ALL());
            const auto d  = subview(A, ip, 1, ALL());
            const auto du = subview(A, ip, 2, ALL());
            thomas(dl, d, du,
                   Kokkos::View<Real*, BulkLayout, Kokkos::MemoryUnmanaged>(
                     X.data() + ip*in.nrow, in.nrow));
          });
        };
        Kokkos::parallel_for(policy, f);
      } else {
        const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
          Kokkos::single(Kokkos::PerTeam(team), [&] () {
            const int ip = team.league_rank();
            const auto dl = subview(A, ip, 0, ALL());
            const auto d  = subview(A, ip, 1, ALL());
            const auto du = subview(A, ip, 2, ALL());
            thomas(dl, d, du, subview(X, ip, ALL(), ALL()));
          });
        };
        Kokkos::parallel_for(policy, f);
      }
      Kokkos::fence();
      t1 = util::gettime();
    }
  } break;
 case Solver::thomas_pack_a1xm: {
   Kokkos::View<Pack***, BulkLayout> Xp("Xp", in.nprob, in.nrow, (in.nrhs + Pack::n - 1)/Pack::n);
   for (int trial = 0; trial < 2; ++trial) {
      Kokkos::deep_copy(A, Acopy);
      pack_data(X, Xp);
      Kokkos::fence();
      t0 = util::gettime();
      const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
        const int ip = team.league_rank();
        const auto dl = subview(A, ip, 0, ALL());
        const auto d  = subview(A, ip, 1, ALL());
        const auto du = subview(A, ip, 2, ALL());
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
          thomas(dl, d, du,
                 Kokkos::View<Pack**, BulkLayout, Kokkos::MemoryUnmanaged>(
                   Xp, ip, ALL(), ALL()));
        });
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      t1 = util::gettime();
    }
    unpack_data(Xp, X);
  } break;
  case Solver::thomas_amxm: {
    // Distinguish between thomas_amxm impls.
    Kokkos::View<Real****> Ap("Ap", in.nprob, 3, in.nrow, in.nrhs);
    for (int trial = 0; trial < 2; ++trial) {
      pack_scalar_matrix(A, Ap, in.nrhs);
      const auto dc = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
          const int i = team.league_rank();
          for (int r = 0; r < in.nrow; ++r)
            for (int c = 0; c < in.nrhs; ++c)
              X(i,r,c) = B(i,r,c);
        });
      };
      Kokkos::parallel_for(policy, dc);
      Kokkos::fence();
      t0 = util::gettime();
      const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
        const int ip = team.league_rank();
        const auto dl = subview(Ap, ip, 0, ALL(), ALL());
        const auto d  = subview(Ap, ip, 1, ALL(), ALL());
        const auto du = subview(Ap, ip, 2, ALL(), ALL());
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
          thomas(dl, d, du, subview(X, ip, ALL(), ALL()));
        });
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      t1 = util::gettime();
    }
  } break;
  case Solver::thomas_pack_amxm: {
    Kokkos::View<Pack***> Xp("Xp", in.nprob, in.nrow, (in.nrhs + Pack::n - 1)/Pack::n);
    Kokkos::View<Pack****> Ap("Ap", in.nprob, 3, in.nrow, (in.nrhs + Pack::n - 1)/Pack::n);
    for (int trial = 0; trial < 2; ++trial) {
      pack_matrix(A, Ap, in.nrhs);
      pack_data(X, Xp);
      Kokkos::fence();
      t0 = util::gettime();
      const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
        const int ip = team.league_rank();
        const auto dl = subview(Ap, ip, 0, ALL(), ALL());
        const auto d  = subview(Ap, ip, 1, ALL(), ALL());
        const auto du = subview(Ap, ip, 2, ALL(), ALL());
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
          thomas(dl, d, du, subview(Xp, ip, ALL(), ALL()));
        });
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      t1 = util::gettime();
    }
    unpack_data(Xp, X);
  } break;
#if defined TRIDIAG_HAVE_LAPACK || defined TRIDIAG_HAVE_MKL
  case Solver::gttr: {
    static const int max_nrow = 512;
    if (in.nrow > max_nrow) throw std::runtime_error("For gttr, run_cpu max nrow is 512.");
    Kokkos::fence();
    for (int trial = 0; trial < 2; ++trial) {
      transpose_data(X, Xt);
      Kokkos::deep_copy(A, Acopy);
      Kokkos::fence();
      t0 = util::gettime();
#     ifdef KOKKOS_ENABLE_OPENMP
#     pragma omp parallel for
#     endif
      for (int ip = 0; ip < in.nprob; ++ip) {
        Real* const p = A.data() + ip*3*in.nrow;
        Real* const dl = p + 1;
        Real* const d  = p +   in.nrow;
        Real* const du = p + 2*in.nrow;
        Real* const x = Xt.data() + ip*in.nrhs*in.nrow;
        int ipiv[max_nrow];
        Real du2[max_nrow];
        gttrf(in.nrow, dl, d, du, du2, ipiv);
        gttrs('n', in.nrow, in.nrhs, dl, d, du, du2, ipiv, x, in.nrow);
      }
      t1 = util::gettime();
    }
    transpose_data(Xt, X);
  } break;
#endif
#ifdef TRIDIAG_HAVE_MKL
  case Solver::dttr: {
    Kokkos::fence();
    for (int trial = 0; trial < 2; ++trial) {
      transpose_data(X, Xt);
      Kokkos::deep_copy(A, Acopy);
      Kokkos::fence();
      t0 = util::gettime();
#     ifdef KOKKOS_ENABLE_OPENMP
#     pragma omp parallel for
#     endif
      for (int ip = 0; ip < in.nprob; ++ip) {
        Real* const p = A.data() + ip*3*in.nrow;
        Real* const dl = p + 1;
        Real* const d  = p +   in.nrow;
        Real* const du = p + 2*in.nrow;
        Real* const x = Xt.data() + ip*in.nrhs*in.nrow;
        dttrfb(in.nrow, dl, d, du);
        dttrsb('n', in.nrow, in.nrhs, dl, d, du, x, in.nrow);
      }
      t1 = util::gettime();
    }
    transpose_data(Xt, X);        
  } break;
#endif
  default:
    std::cout << "run_cpu does not support " << Solver::convert(in.method)
              << "\n";    
  }
  const auto et = t1 - t0;
  printf("run: et %1.3e et/datum %1.3e\n", et, et/(in.nprob*in.nrow*in.nrhs));
  Real re; {
    const auto ip = std::max(0, in.nprob-1);
    const auto dl = subview(Acopy, ip, 0, ALL());
    const auto d  = subview(Acopy, ip, 1, ALL());
    const auto du = subview(Acopy, ip, 2, ALL());
    matvec(dl, d, du,
           subview(X, in.nprob-1, ALL(), ALL()),
           subview(Y, in.nprob-1, ALL(), ALL()));
    re = reldif(
      subview(B, in.nprob-1, ALL(), ALL()),
      subview(Y, in.nprob-1, ALL(), ALL()));
  }
  if (re > 50*std::numeric_limits<Real>::epsilon())
    std::cout << "run: " << " re " << re << "\n";
}

int main (int argc, char** argv) {
  std::cout << "avx:" << util::active_avx_string() << "\n";
#ifdef KOKKOS_ENABLE_OPENMP
  std::cout << "omp: " << omp_get_max_threads() << "\n";
#endif
  int stat = 0;
  Kokkos::initialize(argc, argv); {
    if (argc > 1) {
      Input in;
      stat = in.parse(argc, argv);
      if (stat) {
        if (in.gpu) run_gpu(in);
        else run_cpu(in);
      }
    } else
      test1();  
  } Kokkos::finalize();
  return stat;
}
