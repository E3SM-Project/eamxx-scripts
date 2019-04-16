#include "li_vanilla.hpp"
#include "li_alg.hpp"
#include "li_kokkos.hpp"
#include "li_vect.hpp"
#include "li_common.hpp"

#include <random>
#include <vector>

int main (int argc, char** argv) {

  Kokkos::initialize(argc, argv);

  std::default_random_engine generator;
  std::uniform_int_distribution<int> k_dist(10,100);
  const Real minthresh = 0.000001;
  const int ncol = 10;

  for (int r = 0; r < 1000; ++r) {
    const int km1 = k_dist(generator);
    const int km2 = k_dist(generator);
    vector_2d_t<Real>
      x1(ncol, std::vector<Real>(km1)),
      x2(ncol, std::vector<Real>(km2)),
      y1(ncol, std::vector<Real>(km1)),
      y2_base(ncol, std::vector<Real>(km2)),
      y2_cmp(ncol, std::vector<Real>(km2));

    for (int i = 0; i < ncol; ++i) {
      li::populate_li_input(km1, km2, x1[i].data(), y1[i].data(), x2[i].data());
    }

    typename li::LiKokkos<Real>::template view_2d<Real>
      x1k("x1k", ncol, km1),
      x2k("x2k", ncol, km2),
      y1k("y1k", ncol, km1),
      y2k("y2k", ncol, km2);

    auto x1km = Kokkos::create_mirror_view(x1k);
    auto x2km = Kokkos::create_mirror_view(x2k);
    auto y1km = Kokkos::create_mirror_view(y1k);

    for (int i = 0; i < ncol; ++i) {
      for (int j = 0; j < km1; ++j) {
        x1km(i, j) = x1[i][j];
        y1km(i, j) = y1[i][j];
      }
      for (int j = 0; j < km2; ++j) {
        x2km(i, j) = x2[i][j];
      }
    }

    Kokkos::deep_copy(x1k, x1km);
    Kokkos::deep_copy(x2k, x2km);
    Kokkos::deep_copy(y1k, y1km);

    using LIV = li::LiVect<Real>;
    LIV vect(ncol, km1, km2, minthresh);
    const int km1_pack = vect.km1_pack();
    const int km2_pack = vect.km2_pack();
    using Pack = typename li::LiVect<Real>::Pack;
    typename li::LiVect<Real>::template view_2d<Pack>
      x1kv("x1kv", ncol, km1_pack),
      x2kv("x2kv", ncol, km2_pack),
      y1kv("y1kv", ncol, km1_pack),
      y2kv("y2kv", ncol, km2_pack);

    auto x1kvm = Kokkos::create_mirror_view(x1kv);
    auto x2kvm = Kokkos::create_mirror_view(x2kv);
    auto y1kvm = Kokkos::create_mirror_view(y1kv);

    for (int i = 0; i < ncol; ++i) {
      for (int j = 0; j < km1_pack; ++j) {
        for (int s = 0; s < SCREAM_PACKN; ++s) {
          x1kvm(i, j)[s] = x1[i][j*SCREAM_PACKN + s];
          y1kvm(i, j)[s] = y1[i][j*SCREAM_PACKN + s];
        }
      }
      for (int j = 0; j < km2_pack; ++j) {
        for (int s = 0; s < SCREAM_PACKN; ++s) {
          x2kvm(i, j)[s] = x2[i][j*SCREAM_PACKN + s];
        }
      }
    }

    Kokkos::deep_copy(x1kv, x1kvm);
    Kokkos::deep_copy(x2kv, x2kvm);
    Kokkos::deep_copy(y1kv, y1kvm);

    li::LiVanilla<Real> vanilla(ncol, km1, km2, minthresh);
    li::LiAlg<Real> alg(ncol, km1, km2, minthresh);
    vanilla.setup(x1, x2);
    alg.setup(x1, x2);
    for (int i = 0; i < ncol; ++i) {
      vanilla.lin_interp(x1[i], x2[i], y1[i], y2_base[i], i);
      alg.lin_interp(x1[i], x2[i], y1[i], y2_cmp[i], i);
    }

    using LIK = li::LiKokkos<Real>;
    LIK kokkos(ncol, km1, km2, minthresh);
    kokkos.setup(x1k, x2k);
    Kokkos::parallel_for("lin-interp-ut-kokkos",
                         kokkos.m_policy,
                         KOKKOS_LAMBDA(typename LIK::MemberType const& team_member) {
      const int i = team_member.league_rank();
      kokkos.lin_interp(util::subview(x1k, i),
                        util::subview(x2k, i),
                        util::subview(y1k, i),
                        util::subview(y2k, i),
                        team_member);
    });

    vect.setup(x1kv, x2kv);
    Kokkos::parallel_for("lin-interp-ut-vect",
                         vect.m_policy,
                         KOKKOS_LAMBDA(typename LIV::MemberType const& team_member) {
      const int i = team_member.league_rank();
      vect.lin_interp(util::subview(x1kv, i),
                      util::subview(x2kv, i),
                      util::subview(y1kv, i),
                      util::subview(y2kv, i),
                      team_member);
    });

    auto y2km  = Kokkos::create_mirror_view(y2k);
    auto y2kvm = Kokkos::create_mirror_view(y2kv);
    Kokkos::deep_copy(y2km, y2k);
    Kokkos::deep_copy(y2kvm, y2kv);
    for (int i = 0; i < ncol; ++i) {
      for (int j = 0; j < km2; ++j) {
        micro_require(y2_base[i][j] == y2_cmp[i][j]);
        micro_require(y2_base[i][j] == y2km(i, j));
        micro_require(y2_base[i][j] == y2kvm(i, j / SCREAM_PACKN)[j % SCREAM_PACKN]);
      }
    }
  }

  Kokkos::finalize();

  return 0;
}
