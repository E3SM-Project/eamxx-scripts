#include "types.hpp"
#include "util.hpp"

static Int unittest () {
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

int main (int argc, char** argv) {
  util::initialize();

  Int out = 0;
  Kokkos::initialize(argc, argv); {
    out = unittest();
  } Kokkos::finalize();

  return out;
}
