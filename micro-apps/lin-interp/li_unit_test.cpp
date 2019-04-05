#include "li_vanilla.hpp"
#include "li_alg.hpp"
#include "li_kokkos.hpp"

#include <random>
#include <vector>

int main (int argc, char** argv) {

  Kokkos::initialize(argc, argv);

  std::default_random_engine generator;
  std::uniform_int_distribution<int> k_dist(10,100);
  std::uniform_real_distribution<Real> x_dist(0.0,1.0);
  std::uniform_real_distribution<Real> y_dist(0.0,100.0);
  const Real minthresh = 0.000001;
  const int ncol = 1;

  for (int i = 0; i < 10000; ++i) {
    const int km1 = k_dist(generator);
    const int km2 = k_dist(generator);
    vector_2d_t<Real>
      x1(1, std::vector<Real>(km1)),
      x2(1, std::vector<Real>(km2)),
      y1(1, std::vector<Real>(km1)),
      y2_base(1, std::vector<Real>(km2)),
      y2_cmp(1, std::vector<Real>(km2));

    for (int j = 0; j < km1; ++j) {
      x1[0][j] = x_dist(generator);
      y1[0][j] = y_dist(generator);
    }
    for (int j = 0; j < km2; ++j) {
      x2[0][j] = x_dist(generator);
    }

    std::sort(x1[0].begin(), x1[0].end());
    std::sort(x2[0].begin(), x2[0].end());

    typename li::LiKokkos<Real>::template view_2d<Real>
      x1k("x1k", 1, km1),
      x2k("x2k", 1, km2),
      y1k("y1k", 1, km1),
      y2k("y2k", 1, km2);

    auto x1km = Kokkos::create_mirror_view(x1k);
    auto x2km = Kokkos::create_mirror_view(x2k);
    auto y1km = Kokkos::create_mirror_view(y1k);

    for (int j = 0; j < km1; ++j) {
      x1km(0, j) = x1[0][j];
      y1km(0, j) = y1[0][j];
    }
    for (int j = 0; j < km2; ++j) {
      x2km(0, j) = x2[0][j];
    }

    Kokkos::deep_copy(x1k, x1km);
    Kokkos::deep_copy(x2k, x2km);
    Kokkos::deep_copy(y1k, y1km);

    li::LiVanilla<Real> vanilla(ncol, km1, km2, minthresh);
    vanilla.lin_interp(x1, x2, y1, y2_base);

    li::LiAlg<Real> alg(ncol, km1, km2, minthresh);
    alg.lin_interp(x1, x2, y1, y2_cmp);

    li::LiKokkos<Real> kokkos(ncol, km1, km2, minthresh);
    kokkos.lin_interp(x1k, x2k, y1k, y2k);

    auto y2km = Kokkos::create_mirror_view(y2k);
    for (int j = 0; j < km2; ++j) {
      micro_require(y2_base[0][j] == y2_cmp[0][j]);
      micro_require(y2_base[0][j] == y2km(0, j));
    }
  }

  Kokkos::finalize();

  return 0;
}
