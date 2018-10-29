#include "types.hpp"
#include "util.hpp"
#include "micro_kokkos.hpp"
#include "scream_pack.hpp"

#include <thread>
#include <array>
#include <algorithm>

static Int unittest_team_policy () {
  Int nerr = 0;

  if (util::is_single_precision<double>::value) ++nerr;
  if ( ! util::is_single_precision<float>::value) ++nerr;

  for (int nk: {128, 122, 255, 42}) {
    const int ni = 1000;
    const auto p = util::ExeSpaceUtils<>::get_default_team_policy(ni, nk);
    std::cout << "ni " << ni << " nk " << nk
              << " league size " << p.league_size()
              << " team_size " << p.team_size() << "\n";
    if (p.league_size() != ni) ++nerr;
#ifdef KOKKOS_ENABLE_CUDA
    if (nk == 42) {
      if (p.team_size() != 64) ++nerr;
    } else {
      if (p.team_size() != 128) ++nerr;
    }
#elif defined MIMIC_GPU && defined KOKKOS_ENABLE_OPENMP
    if (omp_get_num_threads() > 1 && p.team_size() == 1) ++nerr;
#endif
  }

  return nerr;
}

static Int unittest_pack () {
  int nerr = 0;
  const int num_bigs = 17;

  using TestBigPack   = scream::pack::Pack<Real, 16>;

  kokkos_1d_t<TestBigPack> test_k_array("test_k_array", num_bigs);
  Kokkos::parallel_reduce("unittest_pack",
                          util::ExeSpaceUtils<>::get_default_team_policy(128, 128),
                          KOKKOS_LAMBDA(const member_type& team, int& total_errs) {

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

namespace unit_test {
struct UnitTest {
static int unittest_workspace()
{
  int nerr = 0;
  const int ints_per_ws = 37;
  static constexpr const int num_ws = 4;
  const int ni = 128;
  const int nk = 128;

  team_policy policy(util::ExeSpaceUtils<>::get_default_team_policy(ni, nk));

  {
    util::WorkspaceManager<double> wsmd(17, num_ws, policy);
    micro_assert(wsmd.m_reserve == 1);
    micro_assert(wsmd.m_size == 17);
  }
  {
    util::WorkspaceManager<char> wsmc(16, num_ws, policy);
    micro_assert(wsmc.m_reserve == 8);
    micro_assert(wsmc.m_size == 16);
    Kokkos::parallel_for(
      "unittest_workspace char", policy,
      KOKKOS_LAMBDA(const member_type& team) {
        auto ws = wsmc.get_workspace(team);
        const auto t1 = ws.take("t1");
        const auto t2 = ws.take("t1");
        ws.release(t1);
        ws.release(t2);
      });
  }
  {
    util::WorkspaceManager<short> wsms(16, num_ws, policy);
    micro_assert(wsms.m_reserve == 4);
    micro_assert(wsms.m_size == 16);
  }

  // Test host-explicit WorkspaceMgr
  {
    team_policy policy_host(util::ExeSpaceUtils<Kokkos::DefaultHostExecutionSpace>::get_default_team_policy(ni, nk));
    util::WorkspaceManager<short, Kokkos::DefaultHostExecutionSpace> wsmh(16, num_ws, policy);
    wsmh.m_data(0, 0) = 0;
  }

  util::WorkspaceManager<int> wsm(ints_per_ws, num_ws, policy);

  micro_assert(wsm.m_reserve == 2);
  micro_assert(wsm.m_size == ints_per_ws);

  Kokkos::parallel_reduce("unittest_workspace", policy, KOKKOS_LAMBDA(const member_type& team, int& total_errs) {

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

      auto ws_dlb = ws.take<double>("doubles");
      if (ws_dlb.extent(0) != 18) ++nerrs_local;
      ws.release(ws_dlb);
    }
    team.team_barrier();

    Kokkos::Array<Unmanaged<kokkos_1d_t<int> >, num_ws> wssub;

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
        Unmanaged<kokkos_1d_t<int> > ws1, ws2, ws3, ws4;
        Kokkos::Array<Unmanaged<kokkos_1d_t<int> >*, num_ws> ptrs = { {&ws1, &ws2, &ws3, &ws4} };
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
          if (wsm.get_index<int>(wssub[w]) != w) ++nerrs_local;
          if (wsm.get_next<int>(wssub[w]) != w+1) ++nerrs_local;
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
              if (!ws.is_active<int>(wssub[w])) ++nerrs_local;
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
};
}

#if 0
static int unittest_team_utils()
{
  int nerr = 0;
  int N = omp_get_max_threads();

  for (int n = 1; n <= N; ++n) {
    const int ni = n*5;
    omp_set_num_threads(n);
    for (int s = 1; s <= n; ++s) {
      const auto p = util::ExeSpaceUtils<>::get_team_policy_force_team_size(ni, s);
      const int c = util::ExeSpaceUtils<>::get_num_concurrent_teams(p);
      util::TeamUtils<> tu(p);
      Kokkos::parallel_reduce("unittest_team_utils", p, KOKKOS_LAMBDA(member_type team_member, int& total_errs) {
        int nerrs_local = 0;
        const int i  = team_member.league_rank();
        int expected_idx = i;
        const int wi = tu.get_workspace_idx(team_member);
        const int thread_num = omp_get_thread_num();

        for (int j = 0; j < n; ++j) {
          if (j == thread_num) {
            if (j == 0) {
              std::cout << "===================================" << std::endl;
            }
            std::cout << "For total_threads: " << n << " and team_size: " << s << ", team: " << i << ", team_rank=" << team_member.team_rank() << ", thread: " << thread_num << " , conc: " << c << ", expected_idx: " << expected_idx << ", idx: " << wi << std::endl;
         }
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
          team_member.team_barrier();
        }

        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team_member, team_member.team_size()), [=] (int t, int& team_errs) {
#ifdef KOKKOS_ENABLE_CUDA
          if (wi != i)  ++team_errs;
#elif defined KOKKOS_ENABLE_OPENMP
          if (wi >= c) ++team_errs;
          if (wi != expected_idx) ++team_errs;
#endif
        }, nerrs_local);
        total_errs += nerrs_local;
        expected_idx += c;
      }, nerr);
      std::cout << "===============================================" << std::endl;
    }
  }

  omp_set_num_threads(N);

  return nerr;
}
#endif

int main (int argc, char** argv) {
  util::initialize();

  int out = 0;
  Kokkos::initialize(argc, argv); {
    out =  unittest_team_policy();
    out += unittest_pack();
    out += unit_test::UnitTest::unittest_workspace();
#if 0
    out += unittest_team_utils();
#endif
  } Kokkos::finalize();

  if (out != 0) std::cout << "Some tests failed" << std::endl;

  return out != 0 ? -1 : 0;
}
