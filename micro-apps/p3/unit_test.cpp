#include "types.hpp"
#include "util.hpp"
#include "micro_kokkos.hpp"
#include "scream_pack.hpp"

#include <thread>
#include <array>
#include <algorithm>

namespace unit_test {

struct UnitWrap {

template <typename D=DefaultDevice>
struct UnitTest {

using MemberType = typename KokkosTypes<D>::MemberType;
using TeamPolicy = typename KokkosTypes<D>::TeamPolicy;
using ExeSpace   = typename KokkosTypes<D>::ExeSpace;

template <typename S>
using view_1d = typename KokkosTypes<D>::template view_1d<S>;
template <typename S>
using view_2d = typename KokkosTypes<D>::template view_2d<S>;
template <typename S>
using view_3d = typename KokkosTypes<D>::template view_3d<S>;

static Int unittest_team_policy () {
  Int nerr = 0;

  if (util::is_single_precision<double>::value) ++nerr;
  if ( ! util::is_single_precision<float>::value) ++nerr;

  for (int nk: {128, 122, 255, 42}) {
    const int ni = 1000;
    const auto p = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ni, nk);
    std::cout << "ni " << ni << " nk " << nk
              << " league size " << p.league_size()
              << " team_size " << p.team_size() << "\n";
    if (p.league_size() != ni) ++nerr;
    if (util::OnGpu<ExeSpace>::value) {
      if (nk == 42) {
        if (p.team_size() != 64) ++nerr;
      } else {
        if (p.team_size() != 128) ++nerr;
      }
    }
    else {
#if defined MIMIC_GPU && defined KOKKOS_ENABLE_OPENMP
      if (omp_get_num_threads() > 1 && p.team_size() == 1) ++nerr;
#endif
    }
  }

  return nerr;
}

static Int unittest_pack () {
  int nerr = 0;
  const int num_bigs = 17;

  using TestBigPack   = scream::pack::Pack<Real, 16>;

  view_1d<TestBigPack> test_k_array("test_k_array", num_bigs);
  Kokkos::parallel_reduce("unittest_pack",
                          util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(128, 128),
                          KOKKOS_LAMBDA(const MemberType& team, int& total_errs) {

    int nerrs_local = 0;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_bigs), [&] (int i) {
      test_k_array(i) = i;
    });

    auto small = scream::pack::repack<4>(test_k_array);
    if (small.extent(0) != 4 * num_bigs) ++nerrs_local;

    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_bigs*4), [&] (int i) {
      for (int p = 0; p < 4; ++p) {
        if (small(i)[p] != i / 4) ++nerrs_local;
      }
    });

    auto big = scream::pack::repack<16>(small);
    if (big.extent(0) != num_bigs) ++nerrs_local;

    total_errs += nerrs_local;
  }, nerr);

  return nerr;
}

static int unittest_workspace()
{
  int nerr = 0;
  const int ints_per_ws = 37;
  static constexpr const int num_ws = 4;
  const int ni = 128;
  const int nk = 128;

  TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ni, nk));

  {
    util::WorkspaceManager<double, D> wsmd(17, num_ws, policy);
    micro_assert(wsmd.m_reserve == 1);
    micro_assert(wsmd.m_size == 17);
  }
  {
    util::WorkspaceManager<char, D> wsmc(16, num_ws, policy);
    micro_assert(wsmc.m_reserve == 8);
    micro_assert(wsmc.m_size == 16);
    Kokkos::parallel_for(
      "unittest_workspace char", policy,
      KOKKOS_LAMBDA(const MemberType& team) {
        auto ws = wsmc.get_workspace(team);
        const auto t1 = ws.take("t1");
        const auto t2 = ws.take("t1");
        ws.release(t1);
        ws.release(t2);
      });
  }
  {
    util::WorkspaceManager<short, D> wsms(16, num_ws, policy);
    micro_assert(wsms.m_reserve == 4);
    micro_assert(wsms.m_size == 16);
  }

  // Test host-explicit WorkspaceMgr
  {
    using HostDevice = Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;
    typename KokkosTypes<HostDevice>::TeamPolicy policy_host(util::ExeSpaceUtils<typename KokkosTypes<HostDevice>::ExeSpace>::get_default_team_policy(ni, nk));
    util::WorkspaceManager<short, HostDevice> wsmh(16, num_ws, policy_host);
    wsmh.m_data(0, 0) = 0; // check on cuda machine
  }

  util::WorkspaceManager<int, D> wsm(ints_per_ws, num_ws, policy);

  micro_assert(wsm.m_reserve == 2);
  micro_assert(wsm.m_size == ints_per_ws);

  Kokkos::parallel_reduce("unittest_workspace", policy, KOKKOS_LAMBDA(const MemberType& team, int& total_errs) {

    int nerrs_local = 0;
    auto ws = wsm.get_workspace(team);

    // Test getting workspaces of different type
    {
      const auto ws_int = ws.take("ints");
      // These nerrs_local increments are write race conditions among threads in
      // a team, but that's OK: nerrs_local doesn't have to be accurate. A 0
      // result will be a true 0 result.
      if (ws_int.extent(0) != ints_per_ws) ++nerrs_local;
      ws.release(ws_int);

      auto ws_dlb = ws.template take<double>("doubles");
      if (ws_dlb.extent(0) != 18) ++nerrs_local;
      ws.release(ws_dlb);
    }
    team.team_barrier();

    Kokkos::Array<Unmanaged<view_1d<int> >, num_ws> wssub;

    // Main test. Test different means of taking and release spaces.
    for (int r = 0; r < 8; ++r) {
      if (r % 4 == 0) {
        for (int w = 0; w < num_ws; ++w) {
          char buf[8] = "ws";
          buf[2] = 48 + w; // 48 is offset to integers in ascii
          wssub[w] = ws.take(buf);
        }
      }
      else {
        Unmanaged<view_1d<int> > ws1, ws2, ws3, ws4;
        Kokkos::Array<Unmanaged<view_1d<int> >*, num_ws> ptrs = { {&ws1, &ws2, &ws3, &ws4} };
        Kokkos::Array<const char*, num_ws> names = { {"ws0", "ws1", "ws2", "ws3"} };
        if (r % 4 == 1) {
          ws.take_many(names, ptrs);
        }
        else if (r % 4 == 2) {
          ws.take_many_contiguous_unsafe(names, ptrs);
        }
        else { // % 4 == 3
          ws.take_many_and_reset(names, ptrs);
        }

        for (int w = 0; w < num_ws; ++w) {
          wssub[w] = *ptrs[w];
        }
      }

      for (int w = 0; w < num_ws; ++w) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ints_per_ws), [&] (Int i) {
          wssub[w](i) = i * w;
        });
      }

      team.team_barrier();

      for (int w = 0; w < num_ws; ++w) {
        // These spaces aren't free, but their metadata should be the same as it
        // was when they were initialized
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
          if (wsm.get_index(wssub[w]) != w) ++nerrs_local;
          if (wsm.get_next(wssub[w]) != w+1) ++nerrs_local;
          char buf[8] = "ws";
          buf[2] = 48 + w; // 48 is offset to integers in ascii
#ifndef NDEBUG
          if (util::strcmp(ws.get_name(wssub[w]), buf) != 0) ++nerrs_local;
          if (ws.get_num_used() != 4) ++nerrs_local;
#endif
          for (int i = 0; i < ints_per_ws; ++i) {
            if (wssub[w](i) != i*w) ++nerrs_local;
          }
        });
      }

      team.team_barrier();

      if (r % 4 == 2) {
        // let take_and_reset do the reset
      }
      else if (r % 2 == 0) {
        ws.reset();
      }
      else {
        for (int w = num_ws - 1; w >= 0; --w) {
          ws.release(wssub[w]);
        }
      }

      team.team_barrier();
    }

#ifndef KOKKOS_ENABLE_CUDA
    // Test weird take/release permutations.
    for (int r = 0; r < 3; ++r) {
      int take_order[]    = {0, 1, 2, 3};
      int release_order[] = {-3, -2, -1, 0};

      do {
        for (int w = 0; w < num_ws; ++w) {
          char buf[8] = "ws";
          buf[2] = 48 + take_order[w]; // 48 is offset to integers in ascii
          wssub[take_order[w]] = ws.take(buf);
        }
        team.team_barrier();

        for (int w = 0; w < num_ws; ++w) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ints_per_ws), [&] (Int i) {
            wssub[w](i) = i * w;
          });
        }

        team.team_barrier();

        // verify stuff
        for (int w = 0; w < num_ws; ++w) {
          Kokkos::single(Kokkos::PerTeam(team), [&] () {
            char buf[8] = "ws";
            buf[2] = 48 + w; // 48 is offset to integers in ascii
#ifndef NDEBUG
            if (util::strcmp(ws.get_name(wssub[w]), buf) != 0) ++nerrs_local;
            if (ws.get_num_used() != 4) ++nerrs_local;
#endif
            for (int i = 0; i < ints_per_ws; ++i) {
              if (wssub[w](i) != i*w) ++nerrs_local;
            }
          });
        }

        team.team_barrier();

        for (int w = 0; w < num_ws; ++w) {
          ws.release(wssub[release_order[w] * -1]);
        }

        team.team_barrier();

        std::next_permutation(release_order, release_order+4);

      } while (std::next_permutation(take_order, take_order+4));
    }
    ws.reset();

    // Test weird take/release permutations.
    for (int r = 0; r < 3; ++r) {
      int actions[] = {-3, -2, -1, 1, 2, 3};
      bool exp_active[] = {false, false, false, false};

      do {
        for (int a = 0; a < 6; ++a) {
          int action = actions[a];
          if (action < 0) {
            action *= -1;
            if (exp_active[action]) {
              ws.release(wssub[action]);
              exp_active[action] = false;
            }
          }
          else {
            if (!exp_active[action]) {
              char buf[8] = "ws";
              buf[2] = 48 + action; // 48 is offset to integers in ascii
              wssub[action] = ws.take(buf);
              exp_active[action] = true;
            }
          }
        }

        for (int w = 0; w < num_ws; ++w) {
          if (exp_active[w]) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ints_per_ws), [&] (Int i) {
              wssub[w](i) = i * w;
            });
          }
        }

        team.team_barrier();

        // verify stuff
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
#ifndef NDEBUG
          int exp_num_active = 0;
#endif
          for (int w = 0; w < num_ws; ++w) {
            char buf[8] = "ws";
            buf[2] = 48 + w; // 48 is offset to integers in ascii
            if (exp_active[w]) {
#ifndef NDEBUG
              if (util::strcmp(ws.get_name(wssub[w]), buf) != 0) ++nerrs_local;
              ++exp_num_active;
              if (!ws.template is_active<int>(wssub[w])) ++nerrs_local;
#endif
              for (int i = 0; i < ints_per_ws; ++i) {
                if (wssub[w](i) != i*w) ++nerrs_local;
              }
            }
          }
#ifndef NDEBUG
          if (ws.get_num_used() != exp_num_active) ++nerrs_local;
#endif
        });

        team.team_barrier();

      } while (std::next_permutation(actions, actions + 6));
    }
    ws.reset();
#endif

    total_errs += nerrs_local;
    team.team_barrier();
  }, nerr);

  wsm.report();

  return nerr;
}

static int unittest_team_utils()
{
  // NOTE: Kokkos does not tolerate changing num_threads post-kokkos-initialization and
  // also does not tolerate multiple initializations in the same process, so this will need
  // to be tested with a bash loop.
  int nerr = 0;
#ifdef KOKKOS_ENABLE_OPENMP
  const int n = omp_get_max_threads();
  const int ni = n*5;
  for (int s = 1; s <= n; ++s) {
    const auto p = util::ExeSpaceUtils<ExeSpace>::get_team_policy_force_team_size(ni, s);
    util::TeamUtils<ExeSpace> tu(p);
    const int c = tu.get_num_concurrent_teams();
    view_2d<int> ws_idxs("ws_idxs", ni, s);
    const int real_ts = omp_get_max_threads() / c;
    std::cout << "thrds " << n << " teamsizeV " << s << " teamsizeR " << real_ts << " ni " << ni << " conc " << c <<  std::endl;
    int kernel_errors = 0;
    Kokkos::parallel_reduce("unittest_team_utils", p, KOKKOS_LAMBDA(MemberType team_member, int& total_errs) {
      int nerrs_local = 0;
      const int i  = team_member.league_rank();
      const int wi = tu.get_workspace_idx(team_member);

#if 0
      const int thread_num = omp_get_thread_num();
      for (int j = 0; j < n; ++j) {
        if (j == thread_num) {
          if (j == 0) {
            std::cout << "===================================" << std::endl;
          }
          std::cout << " For total_threads: " << n << " league size " << team_member.league_size() << " and team_size: " << s << ", team: " << i << ", team_rank=" << team_member.team_rank() << ", thread: " << thread_num << " , conc: " << c << ", idx: " << wi << std::endl;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        team_member.team_barrier();
      }
#endif
      ws_idxs(i, team_member.team_rank()) = wi+1; // never zero
      if (wi >= c) ++nerrs_local;

      total_errs += nerrs_local;
    }, kernel_errors);
#if 0
    std::cout << "===================== DONE ==========================" << std::endl;
#endif
    nerr += kernel_errors;

    // post processing
    const int teams_per_idx = (ni + c - 1) / c;
    for (int i = 0; i < ni; ++i) {
      int exp_wi = 0;
      // all threads in a team should share wsidx
      for (int t = 0; t < s; ++t) {
        int curr_wi = ws_idxs(i, t);
#if 0
        std::cout << "idxs(" << i << ", " << t << ") = " << curr_wi << std::endl;
#endif
        if (t == 0) exp_wi = curr_wi;
        if (curr_wi == 0) ++nerr;
        else if (curr_wi != exp_wi) {
          std::cout << "SHARING ERROR for ws_idxs(" << i << ", " << t << "), " << curr_wi << " != " << exp_wi << std::endl;
          ++nerr;
        }
      }
    }

    // Check that each wsidx used correct number of times
    for (int ci = 1; ci <= c; ++ci) {
      for (int t = 0; t < s; ++t) {
        int cnt = 0;
        for (int i = 0; i < ni; ++i) {
          if (ws_idxs(i,t) == ci) ++cnt;
        }
        if (cnt != teams_per_idx && cnt != teams_per_idx-1) {
          std::cout << "CONC ERROR for ws_idx " << ci << ", was used " << cnt << " times, expected about " << teams_per_idx << "." << std::endl;
          ++nerr;
        }
      }
    }
  }
#endif
  return nerr;
}

};
};
} // namespace unit_test

template <typename D=DefaultDevice>
using UnitTest = unit_test::UnitWrap::UnitTest<D>;

static void expect_another_arg (Int i, Int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}

int main (int argc, char** argv) {
  util::initialize();

  using HostDevice = Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace::memory_space>;

  const auto N =
#ifdef KOKKOS_ENABLE_OPENMP
    omp_get_max_threads() // naturally configurable by machine since env var
#else
    1
#endif
    ;

  int lower = 1, upper=N, increment=1; // defaults
  for (int i = 0; i < argc; ++i) {
    if (util::eq(argv[i], "-l", "--lower")) { expect_another_arg(i, argc); lower     = std::atoi(argv[++i]); }
    if (util::eq(argv[i], "-u", "--upper")) { expect_another_arg(i, argc); upper     = std::atoi(argv[++i]); }
    if (util::eq(argv[i], "-i", "--inc"))   { expect_another_arg(i, argc); increment = std::atoi(argv[++i]); }
  }

  int out = 0; // running error count

  // thread-insensitive tests
  Kokkos::initialize(argc, argv); {
    out += UnitTest<>::unittest_pack();
  } Kokkos::finalize();

  // thread-sensitive tests
  for (int nt = lower; nt <= upper; nt+=increment) { // #threads sweep
#ifdef KOKKOS_ENABLE_OPENMP
    omp_set_num_threads(nt);
#endif
    Kokkos::initialize(argc, argv); {
      out += UnitTest<>::unittest_team_policy();
      out += UnitTest<>::unittest_workspace();
      out += UnitTest<HostDevice>::unittest_team_utils();

#ifdef KOKKOS_ENABLE_CUDA
      // Force host testing on CUDA
      out += UnitTest<HostDevice>::unittest_team_policy();
      out += UnitTest<HostDevice>::unittest_workspace();
#endif
    } Kokkos::finalize();
  }

  if (out != 0) std::cout << "Some tests failed" << std::endl;

  return out != 0 ? -1 : 0;
}
