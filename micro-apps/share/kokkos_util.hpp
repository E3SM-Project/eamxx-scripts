#ifndef KOKKOS_UTIL_HPP
#define KOKKOS_UTIL_HPP

#include "types.hpp"
#include "scream_assert.hpp"
#include "scream_arch.hpp"
#include "micro_kokkos.hpp"

#ifdef _OPENMP
# include <omp.h>
#endif

#include <chrono>

namespace util {

/*
 * Kokkos-related utilities.
 */

/*
 * ExeSpaceUtils is essentially a TeamPolicy factory. TeamPolicy objects
 * are what kokkos uses to define a thread layout (num teams, threads/team)
 * for a parallel kernel. On non-GPU archictures, we will generally have
 * thread teams of 1.
 */
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

/*
 * Specialization of above for Cuda execution space. Many GPU architectures can
 * support a great number of threads, so we'll need to expose additional
 * parallelism by having many threads per team.  This is due to having more
 * threads than the main kernel loop has indices.
 */
#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ExeSpaceUtils<Kokkos::Cuda> {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::Cuda>;

  static TeamPolicy get_default_team_policy (Int ni, Int nk) {
    return TeamPolicy(ni, std::min(128, 32*((nk + 31)/32)));
  }

  static TeamPolicy get_team_policy_force_team_size (Int ni, Int team_size) {
    return TeamPolicy(ni, team_size);
  }
};
#endif

/*
 * TeamUtils contains utilities for getting concurrency info for
 * thread teams. Don't use _TeamUtilsCommonBase directly, use
 * TeamUtils.
 */

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class _TeamUtilsCommonBase
{
 protected:
  int _team_size, _num_teams, _max_threads;

 public:
  template <typename TeamPolicy>
  _TeamUtilsCommonBase(const TeamPolicy& policy) : _team_size(0)
  {
    _max_threads = ExeSpace::concurrency() / ( (!is_single_precision<Real>::value && OnGpu<ExeSpace>::value) ? 2 : 1);
    const int team_size = policy.team_size();
    _num_teams = _max_threads / team_size;
    _team_size = _max_threads / _num_teams;
  }

  // How many thread teams can run concurrently
  int get_num_concurrent_teams() const { return _num_teams; }

  // How many threads can run concurrently
  int get_max_concurrent_threads() const { return _max_threads; }

  // How many ws slots are there
  int get_num_ws_slots() const { return _num_teams; }

  /*
   * Of the C concurrently running teams, which "slot" is open
   * for the given team.
   */
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const
  { return 0; }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  void release_workspace_idx(const MemberType& team_member, int ws_idx) const
  { }
};

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class TeamUtils : public _TeamUtilsCommonBase<ExeSpace>
{
 public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy, const Real& = 1.0) :
    _TeamUtilsCommonBase<ExeSpace>(policy)
  { }
};

/*
 * Specialization for OpenMP execution space
 */
#ifdef KOKKOS_ENABLE_OPENMP
template <>
class TeamUtils<Kokkos::OpenMP> : public _TeamUtilsCommonBase<Kokkos::OpenMP>
{
 public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy, const Real& = 1.0) :
    _TeamUtilsCommonBase<Kokkos::OpenMP>(policy)
  { }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const
  { return omp_get_thread_num() / _team_size; }
};
#endif

/*
 * Specialization for Cuda execution space.
 */
#ifdef KOKKOS_ENABLE_CUDA
template <>
class TeamUtils<Kokkos::Cuda> : public _TeamUtilsCommonBase<Kokkos::Cuda>
{
  using Device = Kokkos::Device<Kokkos::Cuda, typename Kokkos::Cuda::memory_space>;
  using flag_type = int; // this appears to be the smallest type that correctly handles atomic operations
  using view_1d = typename KokkosTypes<Device>::view_1d<flag_type>;
  using RandomGenerator = Kokkos::Random_XorShift64_Pool<Kokkos::Cuda>;
  using rnd_type = typename RandomGenerator::generator_type;

  bool             _need_ws_sharing; // true if there are more teams in the policy than can be run concurrently
  int              _num_ws_slots;    // how many workspace slots (potentially more than the num of concurrent teams due to overprovision factor)
  view_1d         _open_ws_slots;    // indexed by ws-idx, true if in current use, else false
  RandomGenerator _rand_pool;

 public:
  template <typename TeamPolicy>
  TeamUtils(const TeamPolicy& policy, const Real& overprov_factor = 1.25) :
    _TeamUtilsCommonBase<Kokkos::Cuda>(policy),
    _need_ws_sharing(policy.league_size() > _num_teams),
    _num_ws_slots( _need_ws_sharing ? _num_teams * overprov_factor : _num_teams),
    _open_ws_slots("open_ws_slots", _need_ws_sharing ? _num_ws_slots : 0),
    _rand_pool(std::chrono::high_resolution_clock::now().time_since_epoch().count())
  { }

  // How many ws slots are there
  int get_num_ws_slots() const { return _num_ws_slots; }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  int get_workspace_idx(const MemberType& team_member) const
  {
    if (!_need_ws_sharing) {
      return team_member.league_rank();
    }
    else {
      int ws_idx = 0;
      Kokkos::single(Kokkos::PerTeam(team_member), [&] () {
        rnd_type rand_gen = _rand_pool.get_state(team_member.league_rank());
        ws_idx = team_member.league_rank() % _num_ws_slots;
        while (!Kokkos::atomic_compare_exchange_strong(&_open_ws_slots(ws_idx), (flag_type) 0, (flag_type)1)) {
          ws_idx = Kokkos::rand<rnd_type, int>::draw(rand_gen) % _num_ws_slots;
        }
      });

      // broadcast the idx to the team with a simple reduce
      int ws_idx_max_reduce;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, 1), [&] (int, int& ws_idx_max) {
        ws_idx_max = ws_idx;
      }, Kokkos::Max<int>(ws_idx_max_reduce));
      team_member.team_barrier();
      return ws_idx_max_reduce;
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION
  void release_workspace_idx(const MemberType& team_member, int ws_idx) const
  {
    if (_need_ws_sharing) {
      team_member.team_barrier();
      Kokkos::single(Kokkos::PerTeam(team_member), [&] () {
        flag_type volatile* const e = &_open_ws_slots(ws_idx);
        *e = (flag_type)0;
      });
    }
  }
};
#endif

// Get a 1d subview of the i-th dimension of a 2d view
template <typename T, typename ...Parms> KOKKOS_FORCEINLINE_FUNCTION
Unmanaged<Kokkos::View<T*, Parms...> >
subview (const Kokkos::View<T**, Parms...>& v_in, const int i) {
  micro_kassert(v_in.data() != nullptr);
  micro_kassert(i < v_in.extent_int(0));
  micro_kassert(i >= 0);
  return Unmanaged<Kokkos::View<T*, Parms...> >(
    &v_in.impl_map().reference(i, 0), v_in.extent(1));
}

} // namespace util

#endif // SCREAM_KOKKOS_UTILS_HPP
