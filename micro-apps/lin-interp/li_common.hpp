#ifndef LI_COMMON_HPP
#define LI_COMMON_HPP

#include "types.hpp"
#include "util.hpp"
#include "scream_arch.hpp"

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
  Scalar ratio = km1 / static_cast<Scalar>(km2);
  // y is a simple linear function
  for (int k = 0; k < km1; ++k) {
    x1_i[k] = k;
    y1_i[k] = k;
  }
  for (int k = 0; k < km2; ++k) {
    x2_i[k] = k * ratio;
  }
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

  std::vector<Scalar> x1_i(km1), y1_i(km1), x2_i(km2), y2_i(km2);

  populate_li_input(km1, km2, x1_i.data(), y1_i.data(), x2_i.data());

  // This time is thrown out, I just wanted to be able to use auto
  auto start = std::chrono::steady_clock::now();

  for (int r = 0; r < repeat+1; ++r) {

    // re-init
    for (int i = 0; i < ncol; ++i) {
      for (int k = 0; k < km1; ++k) {
        x1[i][k] = x1_i[k];
        y1[i][k] = y1_i[k];
      }
      for (int k = 0; k < km2; ++k) {
        x2[i][k] = x2_i[k];
        y2[i][k] = y2_i[k];
      }
    }

    lik.lin_interp(x1, x2, y1, y2);

    if (r == 0) {
      start = std::chrono::steady_clock::now();
    }
  }

  auto finish = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);

  const double report_time = (1e-6*duration.count()) / repeat;

  printf("Time = %1.3e seconds\n", report_time);

  // Dump the results to a file. This will allow us to do result comparisons between
  // other runs.
  const Scalar* flat_y2 = util::flatten(y2);
  dump_to_file_li(LIK::NAME, flat_y2, ncol, km1, km2, minthresh);
  delete[] flat_y2;
}

} // namespace li

#define common_main(exename)                                            \
  util::initialize();                                                   \
  micro_require_msg(argc == 6, "Usage: " #exename " ncol km1 km2 minthresh repeat"); \
  int ncol(atoi(argv[1])), km1(atoi(argv[2])), km2(atoi(argv[3])), repeat(atoi(argv[5])); \
  Real minthresh(atof(argv[4]))
#endif
