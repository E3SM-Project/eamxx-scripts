#ifndef INCLUDE_UTIL
#define INCLUDE_UTIL

#include "types.hpp"

#include <cstdio>
#include <sstream>
#include <memory>
#include <exception>

#ifndef KOKKOS_ENABLE_CUDA
 #include <cmath>
 #include <algorithm>
#endif

#ifdef _OPENMP
# include <omp.h>
#endif

#ifdef FPE
# include <xmmintrin.h>
#endif

#ifndef NDEBUG
#define micro_assert(condition) do {                                    \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)

#define micro_kernel_assert(condition) do {     \
    if ( ! (condition))                         \
      Kokkos::abort(#condition);                \
  } while (0)

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
#define micro_kernel_throw_if(condition, message) do {              \
    if (condition)                                                  \
      Kokkos::abort(#condition " led to the exception\n" message);  \
  } while (0)


#define common_main(exename)                                                                \
  util::initialize();                                                                       \
  micro_throw_if(argc != 6, "Usage: " #exename "ni nk time_step_len num_steps kdir");       \
  int ni(atoi(argv[1])), nk(atoi(argv[2])), ts(atoi(argv[4])), kdir(atoi(argv[5]));         \
  Real dt(atof(argv[3]));                                                                   \
  micro_throw_if(kdir != -1 && kdir != 1, "kdir must be -1 or 1"); \
  p3::micro_sed::p3_init_cpp<Real>()

namespace util {

struct FILECloser { void operator() (FILE* fh) { fclose(fh); } };
using FILEPtr = std::unique_ptr<FILE, FILECloser>;

template<typename T>
void write (const T* v, size_t sz, const FILEPtr& fid) {
  size_t nwrite = fwrite(v, sizeof(T), sz, fid.get());
  micro_throw_if(nwrite != sz, "write: nwrite = " << nwrite << " sz = " << sz);
}

template<typename T>
void read (T* v, size_t sz, const FILEPtr& fid) {
  size_t nread = fread(v, sizeof(T), sz, fid.get());
  micro_throw_if(nread != sz, "read: nread = " << nread << " sz = " << sz);
}

inline bool
eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

template <typename Real> struct is_single_precision {};
template <> struct is_single_precision<float> { enum : bool { value = true }; };
template <> struct is_single_precision<double> { enum : bool { value = false }; };

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

template <typename Integer> KOKKOS_INLINE_FUNCTION
void set_min_max (const Integer& lim0, const Integer& lim1,
                  Integer& min, Integer& max) {
  min = util::min(lim0, lim1);
  max = util::max(lim0, lim1);
}

inline
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

inline
void dump_arch()
{
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
#ifdef _OPENMP
         omp_get_max_threads()
#else
         1
#endif
         );
}

inline
void initialize()
{
#ifdef FPE
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
                         ~( _MM_MASK_INVALID |
                            _MM_MASK_DIV_ZERO |
                            _MM_MASK_OVERFLOW ));
#endif
}

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
struct ExeSpaceUtils {
  static team_policy get_default_team_policy (Int ni, Int nk) {
#ifdef MIMIC_GPU
    return team_policy(ni, 7);
#else
    return team_policy(ni, 1);
#endif
  }

  template <typename TeamPolicy>
  static int get_num_concurrent_teams(const TeamPolicy& policy)
  {
    const int team_size = policy.team_size();
    const int concurrency = ExeSpace::concurrency();
    return (concurrency + team_size - 1) / team_size;
  }

  template <typename TeamPolicy, typename MemberType>
  static int get_team_idx(const TeamPolicy& policy, const MemberType& team_member)
  {
    return omp_get_thread_num() / policy.team_size();
  }
};

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ExeSpaceUtils<Kokkos::Cuda> {
  static team_policy get_default_team_policy (Int ni, Int nk) {
    return team_policy(ni, std::min(128, 32*((nk + 31)/32)));
  }

  template <typename TeamPolicy>
  static int get_num_concurrent_teams(const TeamPolicy& policy)
  {
    return policy.league_size();
  }

  template <typename TeamPolicy, typename MemberType>
  static int get_team_idx(const TeamPolicy& policy, const MemberType& team_member)
  {
    return team_member.league_rank();
  }
};
#endif

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
struct TeamUtils
{
  int _team_size;

  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy) : _team_size(0)
  {
    const int max_threads = omp_get_max_threads();
    const int team_size = policy.team_size();
    const int num_teams = max_threads / team_size;
    _team_size = max_threads / num_teams;
  }

  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const int team_id) const
  {
    return omp_get_thread_num() / _team_size;
  }
};

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct TeamUtils<Kokkos::Cuda>
{
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy) {}

  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const int team_id) const { return team_id; }
};
#endif

} // namespace util

extern "C" {

void dump_arch_f90();

}

#endif
