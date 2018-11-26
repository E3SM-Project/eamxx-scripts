#ifndef KOKKOS_UTIL_HPP
#define KOKKOS_UTIL_HPP

#include "types.hpp"
#include "scream_assert.hpp"

#ifdef _OPENMP
# include <omp.h>
#endif

namespace util {

/*
 * Kokkos-related utilities.
 *
 * ExeSpaceUtils is essentially a TeamPolicy factory.
 *
 * TeamUtils contains utilities for thread teams
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

#endif
