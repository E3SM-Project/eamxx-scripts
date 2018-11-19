#ifndef INCLUDE_UTIL
#define INCLUDE_UTIL

#include "types.hpp"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <memory>
#include <exception>
#include <map>

#ifndef KOKKOS_ENABLE_CUDA
 #include <cmath>
 #include <algorithm>
#endif

#ifdef _OPENMP
# include <omp.h>
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

#define common_main(exename)                                                                               \
  util::initialize();                                                                                      \
  micro_throw_if(argc != 7, "Usage: " #exename " ni nk time_step_len num_steps kdir repeat");              \
  int ni(atoi(argv[1])), nk(atoi(argv[2])), ts(atoi(argv[4])), kdir(atoi(argv[5])), repeat(atoi(argv[6])); \
  Real dt(atof(argv[3]));                                                                                  \
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

bool eq(const std::string& a, const char* const b1, const char* const b2 = 0);

template <typename Scalar> struct is_single_precision {};
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
KOKKOS_INLINE_FUNCTION
size_t strlen(const char* str)
{
  micro_kernel_assert(str != NULL);
  const char *char_ptr;
  for (char_ptr = str; ; ++char_ptr)  {
    if (*char_ptr == '\0') return char_ptr - str;
  }
}
KOKKOS_INLINE_FUNCTION
void strcpy(char* dst, const char* src)
{
  micro_kernel_assert(dst != NULL && src != NULL);
  while(*dst++ = *src++);
}
KOKKOS_INLINE_FUNCTION
int strcmp(const char* first, const char* second)
{
  while(*first && (*first == *second))
  {
    first++;
    second++;
  }
  return *(const unsigned char*)first - *(const unsigned char*)second;
}
#else
using std::min;
using std::max;
using std::isfinite;
using std::max_element;
using std::strlen;
using std::strcpy;
using std::strcmp;
#endif

template <typename Integer> KOKKOS_INLINE_FUNCTION
void set_min_max (const Integer& lim0, const Integer& lim1,
                  Integer& min, Integer& max) {
  min = util::min(lim0, lim1);
  max = util::max(lim0, lim1);
}

template <typename Integer, typename Integer1> KOKKOS_INLINE_FUNCTION
void set_min_max (const Integer& lim0, const Integer& lim1,
                  Integer& min, Integer& max, const Integer1& vector_size) {
  min = util::min(lim0, lim1) / vector_size;
  max = util::max(lim0, lim1) / vector_size;
}

std::string active_avx_string();

void dump_arch();

void initialize();

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
struct ExeSpaceUtils {
  using TeamPolicy = Kokkos::TeamPolicy<ExeSpace>;

  static TeamPolicy get_default_team_policy (Int ni, Int nk) {
#ifdef MIMIC_GPU
    const int max_threads = ExeSpace::concurrency();
    const int team_size = max_threads < 7 ? max_threads : 7;
    return TeamPolicy(ni, team_size);
#else
    return TeamPolicy(ni, 1);
#endif
  }

  static TeamPolicy get_team_policy_force_team_size (Int ni, Int team_size) {
    return TeamPolicy(ni, team_size);
  }
};

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ExeSpaceUtils<Kokkos::Cuda> {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::Cuda>;

  static TeamPolicy get_default_team_policy (Int ni, Int nk) {
    return TeamPolicy(ni, std::min(128, 32*((nk + 31)/32)));
  }
};
#endif

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class TeamUtils
{
  int _team_size, _num_teams;

public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy) : _team_size(0)
  {
    const int max_threads = ExeSpace::concurrency();
    const int team_size = policy.team_size();
    _num_teams = max_threads / team_size;
    _team_size = max_threads / _num_teams;
  }

  int get_num_concurrent_teams() const { return _num_teams; }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const
  {
    return omp_get_thread_num() / _team_size;
  }
};

#ifdef KOKKOS_ENABLE_CUDA
template <>
class TeamUtils<Kokkos::Cuda>
{
  int _num_teams;

public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy) { _num_teams = policy.league_size(); }

  int get_num_concurrent_teams() const { return _num_teams; }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const { return team_member.league_rank(); }
};
#endif

template <typename Scalar>
void dump_to_file(const char* filename,
                  const Scalar* qr, const Scalar* nr, const Scalar* th, const Scalar* dzq, const Scalar* pres, const Scalar* prt_liq,
                  const int ni, const int nk, const Scalar dt, const int ts, int ldk = -1)
{
  if (ldk < 0) ldk = nk;

  std::string full_fn(filename);
  full_fn += "_perf_run.dat" + std::to_string(sizeof(Scalar));

  FILEPtr fid(fopen(full_fn.c_str(), "w"));
  micro_throw_if( !fid, "dump_to_file can't write " << filename);

  write(&ni, 1, fid);
  write(&nk, 1, fid);
  write(&dt, 1, fid);
  write(&ts, 1, fid);
  // Account for possible alignment padding.
  for (int i = 0; i < ni; ++i) util::write(qr + ldk*i, nk, fid);
  for (int i = 0; i < ni; ++i) util::write(nr + ldk*i, nk, fid);
  for (int i = 0; i < ni; ++i) util::write(th + ldk*i, nk, fid);
  for (int i = 0; i < ni; ++i) util::write(dzq + ldk*i, nk, fid);
  for (int i = 0; i < ni; ++i) util::write(pres + ldk*i, nk, fid);
  write(prt_liq, ni, fid);
}

template <typename ExeSpace>
struct OnGpu { enum : bool { value = false }; };
#ifdef KOKKOS_ENABLE_CUDA
template <> struct OnGpu<Kokkos::Cuda> { enum : bool { value = true }; };
#endif

template <typename T, typename ...Parms> KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<Kokkos::View<T*, Parms...> >
subview (const Kokkos::View<T**, Parms...>& v_in, const int i) {
  micro_kernel_assert(v_in.data() != nullptr);
  micro_kernel_assert(i < v_in.extent_int(0));
  micro_kernel_assert(i >= 0);
  return Unmanaged<Kokkos::View<T*, Parms...> >(
    &v_in.impl_map().reference(i, 0), v_in.extent(1));
}

} // namespace util

extern "C" {

void dump_arch_f90();

}

#endif
