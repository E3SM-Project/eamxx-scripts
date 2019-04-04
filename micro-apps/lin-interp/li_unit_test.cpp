#include "li_vanilla.hpp"
#include "li_alg.hpp"

#include <random>
#include <vector>

int main (int argc, char** argv) {

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

    li::LiVanilla<Real> vanilla(ncol, km1, km2, minthresh);
    vanilla.lin_interp(x1, x2, y1, y2_base);

    li::LiAlg<Real> alg(ncol, km1, km2, minthresh);
    alg.lin_interp(x1, x2, y1, y2_cmp);

    for (int j = 0; j < km2; ++j) {
      micro_require(y2_base[0][j] == y2_cmp[0][j]);
    }
  }

  return 0;
}
