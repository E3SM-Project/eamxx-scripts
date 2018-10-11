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

static int unittest_workspace_1thrd()
{
  int nerr = 0;
  int N = omp_get_max_threads();
  static constexpr const int ints_per_ws = 37;
  static constexpr const int num_ws = 4;

  omp_set_num_threads(1);

  util::WorkSpace ws(sizeof(int) * ints_per_ws, num_ws);

  // Dummy loop just for making a parallel region
  //Kokkos::parallel_for(1, KOKKOS_LAMBDA (const int ignore) {
    std::array<Unmanaged<kokkos_1d_t<int> >, num_ws > wss;

    for (int r = 0; r < 2; ++r) {
      for (int w = 0; w < num_ws; ++w) {
        std::ostringstream oss;
        oss << "ws" << w;
        wss[w] = ws.take<int>(oss.str().c_str());
      }

      for (int w = 0; w < num_ws; ++w) {
        for (int i = 0; i < ints_per_ws; ++i) {
          wss[w](i) = i * w;
        }

        // These spaces aren't free, but their metadata should be the same as it
        // was when they were initialized
        if (ws.get_index<int>(wss[w]) != w) ++nerr;
        if (ws.get_next<int>(wss[w]) != w+1) ++nerr;
      }

      for (int w = num_ws - 1; w >= 0; --w) {
        ws.release<int>(wss[w]);
      }
    }

    wss[0] = ws.take<int>("first");
    wss[1] = ws.take<int>("second");
    wss[2] = ws.take<int>("third");

    ws.release<int>(wss[1]);
    if (ws.get_next<int>(wss[1]) != 3) ++nerr;

    wss[1] = ws.take<int>("second part2");
    if (ws.get_index<int>(wss[1]) != 1) ++nerr;
    //});

  omp_set_num_threads(N);

  return nerr;
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

        for (int j = 0; j < omp_get_num_threads(); ++j) {
          if (j == thread_num) {
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
    out += unittest_workspace_1thrd();
#if 0
    out += unittest_team_utils();
#endif
  } Kokkos::finalize();

  return out;
}
