#include "types.hpp"
#include "util.hpp"
#include "micro_kokkos.hpp"

#include <thread>
#include <array>

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
  util::WorkspaceManager<int> wsm(ints_per_ws, num_ws, policy);

  micro_assert(wsm.m_reserve == 2);
  micro_assert(wsm.m_size == ints_per_ws);

  {
    util::WorkspaceManager<double> wsmd(17, num_ws, policy);
    micro_assert(wsmd.m_reserve == 1);
    micro_assert(wsmd.m_size == 17);
  }
  {
    util::WorkspaceManager<char> wsmc(16, num_ws, policy);
    micro_assert(wsmc.m_reserve == 8);
    micro_assert(wsmc.m_size == 16);
  }
  {
    util::WorkspaceManager<short> wsmc(16, num_ws, policy);
    micro_assert(wsmc.m_reserve == 4);
    micro_assert(wsmc.m_size == 16);
  }

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

    for (int r = 0; r < 2; ++r) {
      {
        for (int w = 0; w < num_ws; ++w) {
          char buf[8] = "ws";
          buf[2] = 48 + w; // 48 is offset to integers in ascii
          wssub[w] = ws.take(buf);
        }
      }

      for (int w = 0; w < num_ws; ++w) {

        char buf[8] = "ws";
        buf[2] = 48 + w; // 48 is offset to integers in ascii
#ifndef NDEBUG
        if (strcmp(ws.get_name(wssub[w]), buf) != 0) ++nerrs_local;
#endif

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ints_per_ws), [&] (Int i) {
          wssub[w](i) = i * w;
        });

        // These spaces aren't free, but their metadata should be the same as it
        // was when they were initialized
        team.team_barrier();
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
          if (wsm.get_index<int>(wssub[w]) != w) ++nerrs_local;
          if (wsm.get_next<int>(wssub[w]) != w+1) ++nerrs_local;

          for (int i = 0; i < ints_per_ws; ++i) {
            if (wssub[w](i) != i*w) ++nerrs_local;
          }
        });
        team.team_barrier();
      }

      {
        for (int w = num_ws - 1; w >= 0; --w) {
          ws.release(wssub[w]);
        }
      }
      team.team_barrier();
    }

    {
      wssub[0] = ws.take("first");
      wssub[1] = ws.take("second");
      wssub[2] = ws.take("third");

      ws.release(wssub[1]);
      if (wsm.get_next<int>(wssub[1]) != 3) ++nerrs_local;

      wssub[1] = ws.take("second part2");
      if (wsm.get_index<int>(wssub[1]) != 1) ++nerrs_local;

      for (int w = num_ws - 2; w >= 0; --w) {
        ws.release(wssub[w]);
      }

      total_errs += nerrs_local;
    }
    team.team_barrier();
  }, nerr);

  // wsm.report(); // uncomment to debug

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

  Int out = 0;
  Kokkos::initialize(argc, argv); {
    out =  unittest_team_policy();
    out += unit_test::UnitTest::unittest_workspace();
#if 0
    out += unittest_team_utils();
#endif
  } Kokkos::finalize();

  if (out != 0) std::cout << "Some tests failed" << std::endl;

  return out;
}
