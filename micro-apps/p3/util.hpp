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

namespace unit_test {
struct UnitTest;
}

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

template <typename Integer, typename Integer1> KOKKOS_INLINE_FUNCTION
void set_min_max (const Integer& lim0, const Integer& lim1,
                  Integer& min, Integer& max, const Integer1& vector_size) {
  min = util::min(lim0, lim1) / vector_size;
  max = util::max(lim0, lim1) / vector_size;
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
    const int max_threads = omp_get_max_threads();
    const int team_size = max_threads < 7 ? max_threads : 7;
    return team_policy(ni, team_size);
#else
    return team_policy(ni, 1);
#endif
  }

  static team_policy get_team_policy_force_team_size (Int ni, Int team_size) {
    return team_policy(ni, team_size);
  }

  template <typename TeamPolicy>
  static int get_num_concurrent_teams(const TeamPolicy& policy)
  {
    const int team_size = policy.team_size();
    const int concurrency = ExeSpace::concurrency();
    return (concurrency + team_size - 1) / team_size;
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

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const
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

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const { return team_member.league_rank(); }
};
#endif

template <typename Real>
void dump_to_file(const char* filename,
                  const Real* qr, const Real* nr, const Real* th, const Real* dzq, const Real* pres, const Real* prt_liq,
                  const int ni, const int nk, const Real dt, const int ts, int ldk = -1)
{
  if (ldk < 0) ldk = nk;

  std::string full_fn(filename);
  full_fn += "_perf_run.dat" + std::to_string(sizeof(Real));

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

template <typename T>
class WorkspaceManager
{
 public:
  WorkspaceManager(int size, int max_used, team_policy policy) :
    m_reserve( (sizeof(T) > 2*sizeof(int)) ? 1 :
               int(std::ceil( (float(sizeof(int))/float(sizeof(T))) * 2 )) ),
    m_size(size + m_reserve),
    m_concurrent_teams(ExeSpaceUtils<>::get_num_concurrent_teams(policy)),
    m_ints_per_ws( (sizeof(T)*m_size) / sizeof(int) ),
    m_tu(policy),
#ifndef NDEBUG
    m_max_used(max_used),
    m_num_used("Workspace.m_num_used", m_concurrent_teams),
#endif
    m_next_slot("Workspace.m_next_slot", m_concurrent_teams),
    m_data("Workspace.m_data", m_concurrent_teams, m_size * max_used)
  {
    // initialize on host
    auto host_mirror = Kokkos::create_mirror_view(m_data);
    int* data = reinterpret_cast<int*>(host_mirror.data());
    for (int t = 0; t < m_concurrent_teams; ++t) {
      for (int i = 0; i < m_max_used; ++i) {
        data[ (m_ints_per_ws * (i+1)) - 2 ] = i; // idx
        data[ (m_ints_per_ws * (i+1)) - 1 ] = i + 1; // next
      }
      data += m_ints_per_ws * m_max_used;
    }

    Kokkos::deep_copy(m_data, host_mirror);
  }

  int get_concurrency() const { return m_concurrent_teams; }

  class Workspace {
   public:
    Workspace(const WorkspaceManager* parent, int ws_idx) : m_parent(parent), m_ws_idx(ws_idx) {}

    template <typename S=T>
    KOKKOS_INLINE_FUNCTION
    Unmanaged<kokkos_1d_t<S> > take(const char* name) const
    { return m_parent->take<S>(name, m_ws_idx); }

    template <typename S=T>
    KOKKOS_INLINE_FUNCTION
    void release(const Unmanaged<kokkos_1d_t<S> >& space) const
    { return m_parent->release<S>(space, m_ws_idx); }

    int index() const { return m_ws_idx; }

   private:
    const WorkspaceManager* m_parent;
    int m_ws_idx;
  };

  Workspace get_workspace(const member_type& team) const
  { return Workspace(this, m_tu.get_workspace_idx(team)); }

 private: // client should be using Workspace

  template <typename S=T>
  KOKKOS_INLINE_FUNCTION
  Unmanaged<kokkos_1d_t<S> > take(const char* name, const int team_idx) const
  {
#ifndef NDEBUG
    micro_kernel_assert(m_num_used(team_idx) < m_max_used);
    m_num_used(team_idx) += 1;
#endif

    auto space = get_space_in_slot<S>(team_idx, m_next_slot(team_idx));
    m_next_slot(team_idx) = get_next<S>(space);

    return space;
  }

  template <typename S=T>
  KOKKOS_INLINE_FUNCTION
  void release(const Unmanaged<kokkos_1d_t<S> >& space, const int team_idx) const
  {
#ifndef NDEBUG
    micro_kernel_assert(m_num_used(team_idx) > 0);
    m_num_used(team_idx) -= 1;
#endif

    set_next<S>(space, m_next_slot(team_idx));
    m_next_slot(team_idx) = get_index<S>(space);
  }

  template <typename S=T>
  KOKKOS_INLINE_FUNCTION
  int get_index(const Unmanaged<kokkos_1d_t<S> >& space) const
  {
    return reinterpret_cast<int*>(space.data())[ m_ints_per_ws - 2 ];
  }

  template <typename S=T>
  KOKKOS_INLINE_FUNCTION
  int get_next(const Unmanaged<kokkos_1d_t<S> >& space) const
  {
    return reinterpret_cast<int*>(space.data())[ m_ints_per_ws - 1 ];
  }

  template <typename S=T>
  KOKKOS_INLINE_FUNCTION
  void set_next(const Unmanaged<kokkos_1d_t<S> >& space, int next) const
  {
    reinterpret_cast<int*>(space.data())[ m_ints_per_ws - 1 ] = next;
  }

  template <typename S=T>
  KOKKOS_INLINE_FUNCTION
  Unmanaged<kokkos_1d_t<S> > get_space_in_slot(const int team_idx, const int slot) const
  {
#ifndef NDEBUG
    micro_kernel_assert(true);
    //int(std::floor(float(m_size - m_reserve) * sizeof(S)/sizeof(T)))
#endif
    if (sizeof(S) == sizeof(T)) {
      return Unmanaged<kokkos_1d_t<T> >( &m_data(team_idx, 0) + m_size*slot, m_size - m_reserve );
    }
    else {
      return Unmanaged<kokkos_1d_t<S> >(
        reinterpret_cast<S*>(&m_data(team_idx, 0) + m_size*slot),
        int((m_size - m_reserve) * float(sizeof(S))/float(sizeof(T)))
                                        );
    }
  }

  friend struct unit_test::UnitTest;

  //
  // data
  //

  int m_reserve, m_size, m_concurrent_teams, m_ints_per_ws;
  util::TeamUtils<> m_tu;
#ifndef NDEBUG
  int m_max_used;
  kokkos_1d_t<int> m_num_used;
#endif
  kokkos_1d_t<int> m_next_slot;
  kokkos_2d_t<T> m_data;
};

} // namespace util

extern "C" {

void dump_arch_f90();

}

#endif
