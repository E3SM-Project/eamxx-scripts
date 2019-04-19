#ifndef LI_COMMON_HPP
#define LI_COMMON_HPP

#include "types.hpp"
#include "util.hpp"
#include "scream_arch.hpp"
#include "kokkos_util.hpp"

#include <random>

extern "C" {

void populate_li_input_from_fortran(int km1, int km2, Real** x1_i, Real** y1_i, Real** x2_i);

bool dump_all_li(const char* filename,
                 const Real** y2,
                 const int ncol, const int km1, const int km2, const Real minthresh);

}

namespace li {

template <typename Scalar>
void populate_li_input(int km1, int km2, Scalar* x1_i, Scalar* y1_i, Scalar* x2_i)
{
  std::default_random_engine generator;
  std::uniform_real_distribution<Real> x_dist(0.0,1.0);
  std::uniform_real_distribution<Real> y_dist(0.0,100.0);

  for (int j = 0; j < km1; ++j) {
    x1_i[j] = x_dist(generator);
    y1_i[j] = y_dist(generator);
  }
  for (int j = 0; j < km2; ++j) {
    x2_i[j] = x_dist(generator);
  }

  // make endpoints same
  x1_i[0] = 0.0;
  x2_i[0] = 0.0;
  x1_i[km1-1] = 1.0;
  x2_i[km2-1] = 1.0;

  std::sort(x1_i, x1_i + km1);
  std::sort(x2_i, x2_i + km2);
}

template <typename Scalar>
void dump_to_file_li(const char* filename, const Scalar* y2, const int ncol, const int km1, const int km2, const Scalar minthresh, int ldk = -1)
{
  if (ldk < 0) ldk = km2;

  std::string full_fn(filename);
  full_fn += "_perf_run.dat" + std::to_string(sizeof(Scalar));

  util::FILEPtr fid(fopen(full_fn.c_str(), "w"));
  micro_require_msg( fid, "dump_to_file can't write " << filename);

  util::write(&ncol, 1, fid);
  util::write(&km1,  1, fid);
  util::write(&km2,  1, fid);
  util::write(&minthresh, 1, fid);

  // Account for possible alignment padding.
  for (int i = 0; i < ncol; ++i) util::write(y2 + ldk*i, km2, fid);
}

template <typename T, typename Scalar, typename D=DefaultDevice>
void populate_kokkos_from_vec(const int num_vert, std::vector<Scalar> const& vec,
                              typename KokkosTypes<D>::template view_1d<T>& device)
{
  const auto mirror = Kokkos::create_mirror_view(device);

  for (int k = 0; k < num_vert; ++k) {
    reinterpret_cast<Scalar*>(mirror.data())[k] = vec[k];
  }

  Kokkos::deep_copy(device, mirror);
}

template <typename Scalar, typename LIK>
void lin_interp_func_wrap(const int ncol, const int km1, const int km2, const Scalar minthresh, const int repeat)
{
  util::dump_arch();

  LIK lik(ncol, km1, km2, minthresh);

  vector_2d_t<Scalar> x1, y1, x2, y2;

  x1.resize(ncol, std::vector<Scalar>(km1));
  y1.resize(ncol, std::vector<Scalar>(km1));
  x2.resize(ncol, std::vector<Scalar>(km2));
  y2.resize(ncol, std::vector<Scalar>(km2));

  std::vector<Scalar> x1_i(km1), y1_i(km1), x2_i(km2);

  populate_li_input(km1, km2, x1_i.data(), y1_i.data(), x2_i.data());

  // init
  for (int i = 0; i < ncol; ++i) {
    for (int k = 0; k < km1; ++k) {
      x1[i][k] = x1_i[k];
      y1[i][k] = y1_i[k];
    }
    for (int k = 0; k < km2; ++k) {
      x2[i][k] = x2_i[k];
      y2[i][k] = 0.0;
    }
  }

  // This time is thrown out, I just wanted to be able to use auto
  auto start = std::chrono::steady_clock::now();

#ifdef LI_TIME_SETUP
  int setup_repeat = repeat;
  int li_repeat = 0;
#else
  int setup_repeat = 0;
  int li_repeat = repeat;
#endif

  for (int r = 0; r < setup_repeat+1; ++r) {

    for (int i = 0; i < ncol; ++i) {
      lik.setup(x1[i], x2[i], i);
    }

#ifdef LI_TIME_SETUP
    if (r == 0) {
      start = std::chrono::steady_clock::now();
    }
#endif
  }

#ifdef LI_TIME_SETUP
  auto finish = std::chrono::steady_clock::now();
#endif

  for (int r = 0; r < li_repeat+1; ++r) {

    for (int i = 0; i < ncol; ++i) {
      lik.lin_interp(x1[i], x2[i], y1[i], y2[i], i);
    }

#ifndef LI_TIME_SETUP
    if (r == 0) {
      start = std::chrono::steady_clock::now();
    }
#endif
  }

#ifndef LI_TIME_SETUP
  auto finish = std::chrono::steady_clock::now();
#endif
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);

  const double report_time = (1e-6*duration.count()) / repeat;
  const char* what =
#ifdef LI_TIME_SETUP
    "setup";
#else
    "li";
#endif

  printf("Time = %1.3e seconds (%s)\n", report_time, what);

  // Dump the results to a file. This will allow us to do result comparisons between
  // other runs.
  const Scalar* flat_y2 = util::flatten(y2);
  dump_to_file_li(LIK::NAME, flat_y2, ncol, km1, km2, minthresh);
  delete[] flat_y2;
}

template <typename Scalar, typename LIK>
void lin_interp_func_wrap_kokkos(const int ncol, const int km1, const int km2, const Scalar minthresh, const int repeat)
{
  util::dump_arch();

  LIK lik(ncol, km1, km2, minthresh);

  using Pack = typename LIK::Pack;

  const int km1_pack = lik.km1_pack();
  const int km2_pack = lik.km2_pack();

  typename LIK::template view_2d<Pack>
    x1("x1", ncol, km1_pack+1),
    y1("y1", ncol, km1_pack+1),
    x2("x2", ncol, km2_pack),
    y2("y2", ncol, km2_pack);

  typename LIK::template view_1d<Pack>
    x1_i("x1_i", km1_pack),
    y1_i("y1_i", km1_pack),
    x2_i("x2_i", km2_pack);

  {
    std::vector<Scalar> x1_iv(km1), y1_iv(km1), x2_iv(km2);
    populate_li_input(km1, km2, x1_iv.data(), y1_iv.data(), x2_iv.data());

    for (auto item : { std::make_pair(&x1_iv, &x1_i), std::make_pair(&y1_iv, &y1_i), std::make_pair(&x2_iv, &x2_i) }) {
      populate_kokkos_from_vec<Pack, Scalar>(item.first->size(), *(item.first), *(item.second));
    }
  }

  int max = std::max(km1_pack, km2_pack);
  Kokkos::parallel_for("init",
                       util::ExeSpaceUtils<typename LIK::ExeSpace>::get_default_team_policy(ncol, max),
                       KOKKOS_LAMBDA(typename LIK::MemberType const& team_member) {
    const int i = team_member.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, max), [=] (int k) {
      if (k < km1_pack) {
        x1(i, k) = x1_i(k);
        y1(i, k) = y1_i(k);
      }
      if (k < km2_pack) {
        x2(i, k) = x2_i(k);
        y2(i, k) = 0.0;
      }
    });
  });

  // This time is thrown out, I just wanted to be able to use auto
  auto start = std::chrono::steady_clock::now();

#ifdef LI_TIME_SETUP
  int setup_repeat = repeat;
  int li_repeat = 0;
#else
  int setup_repeat = 0;
  int li_repeat = repeat;
#endif
  for (int r = 0; r < setup_repeat+1; ++r) {

    Kokkos::parallel_for("setup",
                         lik.m_policy,
                         KOKKOS_LAMBDA(typename LIK::MemberType const& team_member) {
      const int i = team_member.league_rank();
      lik.setup(team_member, util::subview(x1, i), util::subview(x2, i));
    });

#ifdef LI_TIME_SETUP
    if (r == 0) {
      start = std::chrono::steady_clock::now();
    }
#endif
  }

#ifdef LI_TIME_SETUP
  auto finish = std::chrono::steady_clock::now();
#endif

  for (int r = 0; r < li_repeat+1; ++r) {

    Kokkos::parallel_for("lin-interp",
                         lik.m_policy,
                         KOKKOS_LAMBDA(typename LIK::MemberType const& team_member) {
      const int i = team_member.league_rank();
      lik.lin_interp(team_member,
                     util::subview(x1, i),
                     util::subview(x2, i),
                     util::subview(y1, i),
                     util::subview(y2, i));
    });

#ifndef LI_TIME_SETUP
    if (r == 0) {
      start = std::chrono::steady_clock::now();
    }
#endif
  }

#ifndef LI_TIME_SETUP
  auto finish = std::chrono::steady_clock::now();
#endif
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);

  const double report_time = (1e-6*duration.count()) / repeat;
  const char* what =
#ifdef LI_TIME_SETUP
    "setup";
#else
    "li";
#endif

  printf("Time = %1.3e seconds (%s)\n", report_time, what);

  // Dump the results to a file. This will allow us to do result comparisons between
  // other runs.
  const auto mirror_y2 = Kokkos::create_mirror_view(y2);
  Kokkos::deep_copy(mirror_y2, y2);
  const Scalar* flat_y2 = reinterpret_cast<Scalar*>(mirror_y2.data());
  dump_to_file_li(LIK::NAME, flat_y2, ncol, km1, km2, minthresh);
}

} // namespace li

#define common_main(exename)                                            \
  util::initialize();                                                   \
  micro_require_msg(argc == 6, "Usage: " #exename " ncol km1 km2 minthresh repeat"); \
  int ncol(atoi(argv[1])), km1(atoi(argv[2])), km2(atoi(argv[3])), repeat(atoi(argv[5])); \
  Real minthresh(atof(argv[4]))
#endif
