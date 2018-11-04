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

///////////////////////////////////////////////////////////////////////////////
template <typename T, typename D=DefaultDevice>
class WorkspaceManager
///////////////////////////////////////////////////////////////////////////////
{
 public:
  using TeamPolicy = typename KokkosTypes<D>::TeamPolicy;
  using MemberType = typename KokkosTypes<D>::MemberType;
  using ExeSpace   = typename KokkosTypes<D>::ExeSpace;

  template <typename S>
  using view_1d = typename KokkosTypes<D>::template view_1d<S>;
  template <typename S>
  using view_2d = typename KokkosTypes<D>::template view_2d<S>;
  template <typename S>
  using view_3d = typename KokkosTypes<D>::template view_3d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KokkosTypes<D>::template view_1d_ptr_array<S, N>;

  WorkspaceManager(int size, int max_used, TeamPolicy policy) :
    m_tu(policy),
    m_concurrent_teams(m_tu.get_num_concurrent_teams()),
    m_reserve( (sizeof(T) > 2*sizeof(int)) ? 1 :
               (2*sizeof(int) + sizeof(T) - 1)/sizeof(T) ),
    m_size(size),
    m_total(m_size + m_reserve),
    m_max_used(max_used),
#ifndef NDEBUG
    m_num_used("Workspace.m_num_used", m_concurrent_teams),
    m_high_water("Workspace.m_high_water", m_concurrent_teams),
    m_active("Workspace.m_active", m_concurrent_teams, m_max_used),
    m_curr_names("Workspace.m_curr_names", m_concurrent_teams, m_max_used, m_max_name_len),
    m_all_names("Workspace.m_all_names", m_concurrent_teams, m_max_names, m_max_name_len),
    // A name's index in m_all_names is used to index into m_counts
    m_counts("Workspace.m_counts", m_concurrent_teams, m_max_names, 2),
#endif
    m_next_slot("Workspace.m_next_slot", m_pad_factor*m_concurrent_teams),
    m_data(Kokkos::ViewAllocateWithoutInitializing("Workspace.m_data"),
           m_concurrent_teams, m_total * m_max_used)
  {
    init(*this, m_data, m_concurrent_teams, m_max_used, m_total);
  }

  int get_concurrency() const { return m_concurrent_teams; }

  void report() const
  {
#ifndef NDEBUG
    auto host_num_used   = Kokkos::create_mirror_view(m_num_used);
    auto host_high_water = Kokkos::create_mirror_view(m_high_water);
    auto host_all_names  = Kokkos::create_mirror_view(m_all_names);
    auto host_counts     = Kokkos::create_mirror_view(m_counts);

    std::cout << "\nWS usage (capped at " << m_max_used << "): " << std::endl;
    for (int t = 0; t < m_concurrent_teams; ++t) {
      std::cout << "WS " << t << " currently using " << host_num_used(t) << std::endl;
      std::cout << "WS " << t << " high-water " << host_high_water(t) << std::endl;
    }

    std::cout << "\nWS deep analysis" << std::endl;
    std::map<std::string, std::tuple<int, int, int> > ws_usage_map;
    for (int t = 0; t < m_concurrent_teams; ++t) {
      std::cout << "  For wsidx " << t << std::endl;
      for (int n = 0; n < m_max_names; ++n) {
        const char* name = &(host_all_names(t, n, 0));
        if (util::strcmp(name, "") == 0) {
          break;
        }
        else {
          const int takes    = host_counts(t, n, 0);
          const int releases = host_counts(t, n, 1);
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
      std::cout << "Workspace '" << kv.first << "' was used by " << std::get<0>(data) << " wsindices with "
                << std::get<1>(data) << " takes and " << std::get<2>(data) << " releases." << std::endl;
    }
#endif
  }

  class Workspace {
   public:
    KOKKOS_INLINE_FUNCTION
    Workspace(const WorkspaceManager& parent, int ws_idx, const MemberType& team) :
      m_parent(parent), m_team(team), m_ws_idx(ws_idx),
      m_next_slot(parent.m_next_slot(m_pad_factor*ws_idx))
    {}

    template <typename S=T>
    KOKKOS_INLINE_FUNCTION
    Unmanaged<view_1d<S> > take(const char* name) const
    {
#ifndef NDEBUG
      change_num_used(1);
#endif

      const auto space = m_parent.get_space_in_slot<S>(m_ws_idx, m_next_slot);

      // We need a barrier here so get_space_in_slot returns consistent results
      // w/in the team.
      m_team.team_barrier();
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        m_next_slot = m_parent.get_next<S>(space);
#ifndef NDEBUG
        change_indv_meta<S>(space, name);
#endif
      });
      // We need a barrier here so that a subsequent call to take or release
      // starts with the metadata in the correct state.
      m_team.team_barrier();

      return space;
    }

    template <size_t N, typename S=T>
    KOKKOS_INLINE_FUNCTION
    void take_many_contiguous_unsafe(const Kokkos::Array<const char*, N>& names,
                                     const view_1d_ptr_array<S, N>& ptrs) const
    {
#ifndef NDEBUG
      change_num_used(N);
      // Verify contiguous
      for (int n = 0; n < static_cast<int>(N); ++n) {
        const auto space = m_parent.get_space_in_slot<S>(m_ws_idx, m_next_slot + n);
        micro_kernel_assert(m_parent.get_next<S>(space) == m_next_slot + n + 1);
      }
#endif

      for (int n = 0; n < static_cast<int>(N); ++n) {
        const auto space = m_parent.get_space_in_slot<S>(m_ws_idx, m_next_slot+n);
        *ptrs[n] = space;
      }

      // We need a barrier here so get_space_in_slot above returns consistent results
      // w/in the team.
      m_team.team_barrier();
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        m_next_slot += N;
#ifndef NDEBUG
        for (int n = 0; n < static_cast<int>(N); ++n) {
          change_indv_meta<S>(*ptrs[n], names[n]);
        }
#endif
      });
      // We need a barrier here so that a subsequent call to take or release
      // starts with the metadata in the correct state.
      m_team.team_barrier();
    }

    template <size_t N, typename S=T>
    KOKKOS_INLINE_FUNCTION
    void take_many(const Kokkos::Array<const char*, N>& names,
                   const view_1d_ptr_array<S, N>& ptrs) const
    {
#ifndef NDEBUG
      change_num_used(N);
#endif

      int next_slot = m_next_slot;
      for (int n = 0; n < static_cast<int>(N); ++n) {
        auto& space = *ptrs[n];
        space = m_parent.get_space_in_slot<S>(m_ws_idx, next_slot);
        next_slot = m_parent.get_next<S>(space);
      }

      // We need a barrier here so get_space_in_slot above returns consistent results
      // w/in the team.
      m_team.team_barrier();
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        m_next_slot = next_slot;
#ifndef NDEBUG
        for (int n = 0; n < static_cast<int>(N); ++n) {
          change_indv_meta<S>(*ptrs[n], names[n]);
        }
#endif
      });
      // We need a barrier here so that a subsequent call to take or release
      // starts with the metadata in the correct state.
      m_team.team_barrier();
    }

    template <size_t N, typename S=T>
    KOKKOS_INLINE_FUNCTION
    void take_many_and_reset(const Kokkos::Array<const char*, N>& names,
                             const view_1d_ptr_array<S, N>& ptrs) const
    {
#ifndef NDEBUG
      change_num_used(N - m_parent.m_num_used(m_ws_idx));
#endif

      for (int n = 0; n < static_cast<int>(N); ++n) {
        const auto space = m_parent.get_space_in_slot<S>(m_ws_idx, n);
        *ptrs[n] = space;
      }

      // We only need to reset the metadata for spaces that are being left free
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(m_team, m_parent.m_max_used - N), [&] (int i) {
          m_parent.init_metadata(m_ws_idx, i+N);
        });

      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        m_next_slot = N;
#ifndef NDEBUG
        // Mark all old spaces as released
        for (int a = 0; a < m_parent.m_max_used; ++a) {
          if (m_parent.m_active(m_ws_idx, a)) {
            change_indv_meta<S>(m_parent.get_space_in_slot<S>(m_ws_idx, a), "", true);
          }
        }

        // Mark all new spaces as taken
        for (int n = 0; n < static_cast<int>(N); ++n) {
          change_indv_meta<S>(*ptrs[n], names[n]);
        }
#endif
      });
      // We need a barrier here so that a subsequent call to take or release
      // starts with the metadata in the correct state.
      m_team.team_barrier();
    }

    // Wrapper so caller doesn't have to specify scalar type.
    template <typename View>
    KOKKOS_FORCEINLINE_FUNCTION
    void release(const View& space, std::enable_if<View::rank == 1>* = 0) const
    { release_impl<typename View::value_type>(space); }

#ifndef NDEBUG
    template <typename View>
    KOKKOS_INLINE_FUNCTION
    const char* get_name(const View& space, std::enable_if<View::rank == 1>* = 0) const
    { return get_name_impl<typename View::value_type>(space); }
#endif

    int index() const { return m_ws_idx; }

    KOKKOS_INLINE_FUNCTION
    void reset() const
    {
      m_team.team_barrier();
#ifndef NDEBUG
      change_num_used(-m_parent.m_num_used(m_ws_idx));
#endif
      m_next_slot = 0;
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(m_team, m_parent.m_max_used), [&] (int i) {
          m_parent.init_metadata(m_ws_idx, i);
        });

#ifndef NDEBUG
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        // Mark all old spaces as released
        for (int a = 0; a < m_parent.m_max_used; ++a) {
          if (m_parent.m_active(m_ws_idx, a)) {
            change_indv_meta<T>(m_parent.get_space_in_slot<T>(m_ws_idx, a), "", true);
          }
        }
      });
#endif

      m_team.team_barrier();
    }

    // Print the linked list. Obviously not a device function.
    void print () {
      m_team.team_barrier();
      Kokkos::single(
        Kokkos::PerTeam(m_team), [&] () {
          std::stringstream ss;
          ss << m_ws_idx << ":";
          auto space = m_parent.get_space_in_slot<T>(m_ws_idx, m_next_slot);
          for (int cnt = 0, nmax = m_parent.m_max_used;
               cnt < nmax;
               ++cnt) {
            ss << " (" << m_parent.get_index<T>(space) << ", "
               << m_parent.get_next<T>(space) << ")";
            space = m_parent.get_space_in_slot<T>(m_ws_idx, m_parent.get_next<T>(space));
          }
          ss << "\n";
          std::cout << ss.str();
        });
    }

   private:
    const WorkspaceManager& m_parent;
    const MemberType& m_team;
    const int m_ws_idx;
    int& m_next_slot;

#ifndef NDEBUG
    template <typename S>
    KOKKOS_INLINE_FUNCTION
    const char* get_name_impl(const Unmanaged<view_1d<S> >& space) const
    {
      const int slot = m_parent.get_index<S>(space);
      return &(m_parent.m_curr_names(m_ws_idx, slot, 0));
    }

    KOKKOS_INLINE_FUNCTION
    void change_num_used(int change_by) const
    {
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        int curr_used = m_parent.m_num_used(m_ws_idx) += change_by;
        micro_kernel_assert(curr_used <= m_parent.m_max_used);
        micro_kernel_assert(curr_used >= 0);
        if (curr_used > m_parent.m_high_water(m_ws_idx)) {
          m_parent.m_high_water(m_ws_idx) = curr_used;
        }
      });
    }

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    void change_indv_meta(const Unmanaged<view_1d<S> >& space, const char* name, bool release=false) const
    {
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        const int slot = m_parent.get_index<S>(space);
        if (!release) {
          micro_kernel_assert(util::strlen(name) < m_max_name_len); // leave one char for null terminator
          micro_kernel_assert(util::strlen(name) > 0);
          micro_kernel_assert(!m_parent.m_active(m_ws_idx, slot));
          char* val = &(m_parent.m_curr_names(m_ws_idx, slot, 0));
          util::strcpy(val, name);
        }
        else {
          micro_kernel_assert(m_parent.m_active(m_ws_idx, slot));
          name = get_name(space);
        }
        const int name_idx = get_name_idx(name, !release);
        const int count_idx = release ? 1 : 0;
        m_parent.m_counts(m_ws_idx, name_idx, count_idx) += 1;
        m_parent.m_active(m_ws_idx, slot) = !release;
      });
    }

    KOKKOS_INLINE_FUNCTION
    int get_alloc_count(const char* name) const
    { return m_parent.m_counts(m_ws_idx, get_name_idx(name), 0); }

    KOKKOS_INLINE_FUNCTION
    int get_release_count(const char* name) const
    { return m_parent.m_counts(m_ws_idx, get_name_idx(name), 1); }

    KOKKOS_INLINE_FUNCTION
    int get_num_used() const
    { return m_parent.m_num_used(m_ws_idx); }

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    bool is_active(const Unmanaged<view_1d<S> >& space) const
    { return m_parent.m_active(m_ws_idx, m_parent.get_index<S>(space));}

    KOKKOS_INLINE_FUNCTION
    int get_name_idx(const char* name, bool add=false) const
    {
      int name_idx = -1;
      for (int n = 0; n < m_max_names; ++n) {
        char* old_name = &(m_parent.m_all_names(m_ws_idx, n, 0));
        if (util::strcmp(old_name, name) == 0) {
          name_idx = n;
          break;
        }
        else if (add && util::strcmp(old_name, "") == 0) {
          util::strcpy(old_name, name);
          name_idx = n;
          break;
        }
      }
      micro_kernel_assert(name_idx != -1);
      return name_idx;
    }
#endif

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    void release_impl(const Unmanaged<view_1d<S> >& space) const
    {
#ifndef NDEBUG
      change_num_used(-1);
      change_indv_meta<S>(space, "", true);
#endif

      // We don't need a barrier before this block b/c it's OK for metadata to
      // change while some threads in the team are still using the bulk data.
      Kokkos::single(Kokkos::PerTeam(m_team), [&] () {
        m_next_slot = m_parent.set_next_and_get_index<S>(space, m_next_slot);
      });
      m_team.team_barrier();
    }

    friend struct unit_test::UnitTest;
  }; // class Workspace

  KOKKOS_INLINE_FUNCTION
  Workspace get_workspace(const MemberType& team) const
  { return Workspace(*this, m_tu.get_workspace_idx(team), team); }

 public: // for Cuda

  static void init(const WorkspaceManager& wm, const view_2d<T>& data,
                   const int concurrent_teams, const int max_used, const int total)
  {
    Kokkos::parallel_for(
      "WorkspaceManager ctor",
      util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(concurrent_teams, max_used),
      KOKKOS_LAMBDA(const MemberType& team) {
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, max_used), [&] (int i) {
            wm.init_metadata(team.league_rank(), i);
          });
      });
  }

 private: // client should be using Workspace

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int get_index(const Unmanaged<view_1d<S> >& space) const
  {
    return reinterpret_cast<const int*>(reinterpret_cast<const T*>(space.data()) - m_reserve)[0];
  }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int get_next(const Unmanaged<view_1d<S> >& space) const
  {
    return reinterpret_cast<const int*>(reinterpret_cast<const T*>(space.data()) - m_reserve)[1];
  }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int set_next_and_get_index(const Unmanaged<view_1d<S> >& space, int next) const
  {
    const auto metadata = reinterpret_cast<int*>(reinterpret_cast<T*>(space.data()) - m_reserve);
    metadata[1] = next;
    return metadata[0];
  }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  Unmanaged<view_1d<S> > get_space_in_slot(const int team_idx, const int slot) const
  {
    return Unmanaged<view_1d<S> >(
      reinterpret_cast<S*>(&m_data(team_idx, slot*m_total) + m_reserve),
      sizeof(T) == sizeof(S) ?
      m_size :
      (m_size*sizeof(T))/sizeof(S));
  }

  KOKKOS_INLINE_FUNCTION
  void init_metadata(const int ws_idx, const int slot) const
  {
    int* const metadata = reinterpret_cast<int*>(&m_data(ws_idx, slot*m_total));
    metadata[0] = slot;     // idx
    metadata[1] = slot + 1; // next
  }

  friend struct unit_test::UnitTest;

  //
  // data
  //

  enum { m_pad_factor =
#ifdef KOKKOS_ENABLE_CUDA  // TODO: Replace with OnGpu<ES>::value
         1,
#else
         32,
#endif
         m_max_name_len = 128,
         m_max_names = 256
  };

  util::TeamUtils<ExeSpace> m_tu;
  int m_concurrent_teams, m_reserve, m_size, m_total, m_max_used;
#ifndef NDEBUG
  view_1d<int> m_num_used;
  view_1d<int> m_high_water;
  view_2d<bool> m_active;
  view_3d<char> m_curr_names;
  view_3d<char> m_all_names;
  view_3d<int> m_counts;
#endif
  view_1d<int> m_next_slot;
  view_2d<T> m_data;
};

} // namespace util

extern "C" {

void dump_arch_f90();

}

#endif
