#ifndef INCLUDE_UTIL
#define INCLUDE_UTIL

#include "types.hpp"

#include <cstdio>
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

///////////////////////////////////////////////////////////////////////////////
template <typename T>
class WorkspaceManager
///////////////////////////////////////////////////////////////////////////////
{
 public:
  WorkspaceManager(int size, int max_used, team_policy policy) :
    m_reserve( (sizeof(T) > 2*sizeof(int)) ? 1 :
               (2*sizeof(int) + sizeof(T) - 1)/sizeof(T) ),
    m_size(size),
    m_total(m_size + m_reserve),
    m_concurrent_teams(ExeSpaceUtils<>::get_num_concurrent_teams(policy)),
    m_max_used(max_used + 4), // a little padding between threads' data
    m_tu(policy),
#ifndef NDEBUG
    m_num_used("Workspace.m_num_used", m_concurrent_teams),
    m_high_water("Workspace.m_high_water", m_concurrent_teams),
 #ifndef KOKKOS_ENABLE_CUDA
    m_curr_names("Workspace.m_curr_names", m_concurrent_teams, m_max_used, 128), // 128 is max name len
    m_all_names("Workspace.m_all_names", policy.league_size(), 1000, 128), // up to 1000 unique names
    m_counts("Workspace.m_counts", policy.league_size(), 1000),
 #endif
#endif
    m_next_slot("Workspace.m_next_slot", m_pad_factor*m_concurrent_teams),
    m_data(Kokkos::ViewAllocateWithoutInitializing("Workspace.m_data"),
           m_concurrent_teams, m_total * m_max_used)
  {
    init(m_data, m_concurrent_teams, m_max_used, m_total);
  }

  int get_concurrency() const { return m_concurrent_teams; }

  void report() const
  {
#ifndef NDEBUG
    auto host_num_used   = Kokkos::create_mirror_view(m_num_used);
    auto host_high_water = Kokkos::create_mirror_view(m_high_water);
 #ifndef KOKKOS_ENABLE_CUDA
    auto host_all_names  = Kokkos::create_mirror_view(m_all_names);
    auto host_counts     = Kokkos::create_mirror_view(m_counts);
 #endif

    std::cout << "\nWS usage (capped at " << m_max_used << "): " << std::endl;
    for (int t = 0; t < m_concurrent_teams; ++t) {
      std::cout << "WS " << t << " currently using " << host_num_used(t) << std::endl;
      std::cout << "WS " << t << " high-water " << host_high_water(t) << std::endl;
    }

 #ifndef KOKKOS_ENABLE_CUDA
    std::cout << "\nWS deep analysis" << std::endl;
    std::map<std::string, std::tuple<int, int, int> > ws_usage_map;
    const int league_size = m_all_names.extent(0);
    for (int t = 0; t < league_size; ++t) {
      std::cout << "  For team " << t << std::endl;
      for (int n = 0; n < 1000; ++n) {
        const char* name = &(host_all_names(t, n, 0));
        if (strcmp(name, "") == 0) {
          break;
        }
        else {
          const int takes    = host_counts(t, n).first;
          const int releases = host_counts(t, n).second;
          std::cout << "    workspace '" << name << "' was taken " << takes
                    << " times and released " << releases << " times" << std::endl;
          if (takes != releases) {
            std::cout << "      POSSIBLE LEAK" << std::endl;
          }
          std::string sname(name);
          if (ws_usage_map.find(sname) == ws_usage_map.end()) {
            ws_usage_map[sname] = std::make_tuple(1, takes, releases);
          }
          else {
            std::get<0>(ws_usage_map[sname]) += 1;
            std::get<1>(ws_usage_map[sname]) += takes;
            std::get<2>(ws_usage_map[sname]) += releases;
          }
        }
      }
    }

    std::cout << "\nWS workspace summary" << std::endl;
    for (auto& kv : ws_usage_map) {
      auto data = kv.second;
      std::cout << "Workspace '" << kv.first << "' was used by " << std::get<0>(data) << " teams with "
                << std::get<1>(data) << " takes and " << std::get<2>(data) << " releases." << std::endl;
    }
 #endif
#endif
  }

  class Workspace {
   public:
    KOKKOS_INLINE_FUNCTION
    Workspace(const WorkspaceManager& parent, int ws_idx, const member_type& team) :
      m_parent(parent), m_team(team), m_ws_idx(ws_idx),
      m_next_slot(parent.m_next_slot(m_pad_factor*ws_idx))
    {}

    template <typename S=T>
    KOKKOS_INLINE_FUNCTION
    Unmanaged<kokkos_1d_t<S> > take(const char* name) const
    {
#ifndef NDEBUG
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        micro_kernel_assert(m_parent.m_num_used(m_ws_idx) < m_parent.m_max_used);
        m_parent.m_num_used(m_ws_idx) += 1;
        if (m_parent.m_num_used(m_ws_idx) > m_parent.m_high_water(m_ws_idx)) {
          m_parent.m_high_water(m_ws_idx) = m_parent.m_num_used(m_ws_idx);
        }
      });
#endif

      const auto space = m_parent.get_space_in_slot<S>(m_ws_idx, m_next_slot);

      // We need a barrier here so get_space_in_slot returns consistent results
      // w/in the team.
      m_team.team_barrier();
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        m_next_slot = m_parent.get_next<S>(space);
#ifndef NDEBUG
 #ifndef KOKKOS_ENABLE_CUDA
        set_name<S>(space, name);
        const int team_rank = m_team.league_rank();
        int name_idx = -1;
        for (int n = 0; n < 1000; ++n) {
          char* old_name = &(m_parent.m_all_names(team_rank, n, 0));
          if (strcmp(old_name, name) == 0) {
            name_idx = n;
            break;
          }
          else if (strcmp(old_name, "") == 0) {
            strcpy(old_name, name);
            name_idx = n;
            break;
          }
        }
        micro_assert(name_idx != -1);
        m_parent.m_counts(team_rank, name_idx).first += 1;
 #endif
#endif
      });
      // We need a barrier here so that a subsequent call to take or release
      // starts with the metadata in the correct state.
      m_team.team_barrier();

      return space;
    }

    // Wrapper so caller doesn't have to specify scalar type.
    template <typename View>
    KOKKOS_FORCEINLINE_FUNCTION
    void release(const View& space, std::enable_if<View::rank == 1>* = 0) const
    { release_impl<typename View::value_type>(space); }

#ifndef NDEBUG
 #ifndef KOKKOS_ENABLE_CUDA
    template <typename View>
    KOKKOS_INLINE_FUNCTION
    const char* get_name(const View& space, std::enable_if<View::rank == 1>* = 0) const
    { return get_name_impl<typename View::value_type>(space); }
 #endif
#endif

    int index() const { return m_ws_idx; }

    KOKKOS_INLINE_FUNCTION
    void reset() const
    {
      m_team.team_barrier();
      m_next_slot = 0;
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(m_team, m_parent.m_max_used), [&] (int i) {
          int* const metadata = reinterpret_cast<int*>(&m_parent.m_data(m_ws_idx, i*m_parent.m_total));
          metadata[0] = i;     // idx
          metadata[1] = i + 1; // next
        });
      m_team.team_barrier();
    }

   private:
    const WorkspaceManager& m_parent;
    const member_type& m_team;
    const int m_ws_idx;
    int& m_next_slot;

#ifndef NDEBUG
 #ifndef KOKKOS_ENABLE_CUDA
    template <typename S>
    KOKKOS_INLINE_FUNCTION
    const char* get_name_impl(const Unmanaged<kokkos_1d_t<S> >& space) const
    {
      const int slot = m_parent.get_index<S>(space);
      return &(m_parent.m_curr_names(m_ws_idx, slot, 0));
    }

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    void set_name(const Unmanaged<kokkos_1d_t<S> >& space, const char* name) const
    {
      micro_kernel_assert(strlen(name) < 128); // leave one char for null terminator
      const int slot = m_parent.get_index<S>(space);
      char* val = &(m_parent.m_curr_names(m_ws_idx, slot, 0));
      strcpy(val, name);
    }
 #endif
#endif

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    void release_impl(const Unmanaged<kokkos_1d_t<S> >& space) const
    {
#ifndef NDEBUG
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        micro_kernel_assert(m_parent.m_num_used(m_ws_idx) > 0);
        m_parent.m_num_used(m_ws_idx) -= 1;
 #ifndef KOKKOS_ENABLE_CUDA
        const char* name = get_name(space);
        int name_idx = -1;
        const int team_rank = m_team.league_rank();
        for (int n = 0; n < 1000; ++n) {
          const char* old_name = &(m_parent.m_all_names(team_rank, n, 0));
          if (strcmp(old_name, name) == 0) {
            name_idx = n;
            break;
          }
        }
        micro_assert(name_idx != -1);
        m_parent.m_counts(team_rank, name_idx).second += 1;
 #endif
      });
#endif

      // We don't need a barrier before this block b/c it's OK for metadata to
      // change while some threads in the team are still using the bulk data.
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        m_next_slot = m_parent.set_next_and_get_index<S>(space, m_next_slot);
      });
      m_team.team_barrier();
    }
  };

  KOKKOS_INLINE_FUNCTION
  Workspace get_workspace(const member_type& team) const
  { return Workspace(*this, m_tu.get_workspace_idx(team), team); }


 public: // for Cuda

  static void init (const kokkos_2d_t<T>& data, const int concurrent_teams,
                    const int max_used, const int total)
  {
    Kokkos::parallel_for(
      "WorkspaceManager ctor",
      util::ExeSpaceUtils<>::get_default_team_policy(concurrent_teams, max_used),
      KOKKOS_LAMBDA(const member_type& team) {
        const int t = team.league_rank();
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, max_used), [&] (int i) {
            int* const metadata = reinterpret_cast<int*>(&data(t, i*total));
            metadata[0] = i;     // idx
            metadata[1] = i + 1; // next
          });
      });
  }

 private: // client should be using Workspace

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int get_index(const Unmanaged<kokkos_1d_t<S> >& space) const
  {
    return reinterpret_cast<const int*>(reinterpret_cast<const T*>(space.data()) - m_reserve)[0];
  }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int get_next(const Unmanaged<kokkos_1d_t<S> >& space) const
  {
    return reinterpret_cast<const int*>(reinterpret_cast<const T*>(space.data()) - m_reserve)[1];
  }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int set_next_and_get_index(const Unmanaged<kokkos_1d_t<S> >& space, int next) const
  {
    const auto metadata = reinterpret_cast<int*>(reinterpret_cast<T*>(space.data()) - m_reserve);
    metadata[1] = next;
    return metadata[0];
  }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  Unmanaged<kokkos_1d_t<S> > get_space_in_slot(const int team_idx, const int slot) const
  {
    return Unmanaged<kokkos_1d_t<S> >(
      reinterpret_cast<S*>(&m_data(team_idx, slot*m_total) + m_reserve),
      sizeof(T) == sizeof(S) ?
      m_size :
      (m_size*sizeof(T))/sizeof(S));
  }

  friend struct unit_test::UnitTest;

  //
  // data
  //

  enum { m_pad_factor =
#ifdef KOKKOS_ENABLE_CUDA  // TODO: Replace with OnGpu<ES>::value
         1
#else
         32
#endif
  };

  int m_reserve, m_size, m_total, m_concurrent_teams, m_max_used;
  util::TeamUtils<> m_tu;
#ifndef NDEBUG
  kokkos_1d_t<int> m_num_used;
  kokkos_1d_t<int> m_high_water;
 #ifndef KOKKOS_ENABLE_CUDA
  kokkos_3d_t<char> m_curr_names; // TODO None of this runs on GPU, so just use std::vector<std::string> type stuff.
  kokkos_3d_t<char> m_all_names;
  kokkos_2d_t<std::pair<int, int> > m_counts;
 #endif
#endif
  kokkos_1d_t<int> m_next_slot;
  kokkos_2d_t<T> m_data;
};

} // namespace util

extern "C" {

void dump_arch_f90();

}

#endif
