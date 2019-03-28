#include "types.hpp"
#include "util.hpp"
#include "scream_arch.hpp"
#include "../micro-sed/cmp.hpp"

#include <vector>
#include <iostream>
#include <exception>

/*
 * A small command-line utility for comparing two sets of results (in the form
 * of a binary file) from two li runs. This is used in our ctests to make
 * sure all the implementations are producing similar results.
 */

using namespace cmp;

int compare_files(const std::string& file1_fn, const std::string& file2_fn, Real tol, bool verbose)
{

  util::FILEPtr fid1(fopen(file1_fn.c_str(), "r")), fid2(fopen(file2_fn.c_str(), "r"));

  int ncol = compare_property<int>(file1_fn, file2_fn, fid1, fid2, "ncol", verbose);
  compare_property<int>(file1_fn, file2_fn, fid1, fid2, "km1", verbose);
  int km2  = compare_property<int>(file1_fn, file2_fn, fid1, fid2, "km2", verbose);
  compare_property<Real>(file1_fn, file2_fn, fid1, fid2, "minthresh", verbose);

  const int size = ncol * km2;
  std::vector<Real> vec1(size), vec2(size);
  util::read(vec1.data(), size, fid1);
  util::read(vec2.data(), size, fid2);
  int bad_cnt = cmp::compare("y2", vec1.data(), vec2.data(), size, tol, verbose);
  if (bad_cnt == 0 && verbose) std::cout << "y2 matched " << size << " items." << std::endl;

  return bad_cnt;
}

int main (int argc, char** argv) {
  util::initialize();

  if (argc < 3) {
    std::cout <<
      argv[0] << " [options] file1 file2\n"
      "Options:\n"
      "  -v        Run with extra verbose output.\n"
      "  -t <tol>  Tolerance for relative error.\n";
    return -1;
  }

  // When performance testing, we're using optimized flags, etc., so
  // even in double precision we have to accept *some* diffs.
  Real tol = util::is_single_precision<Real>::value ? 2e-5 : 1e-14;
  bool verbose = false;
  for (Int i = 1; i < argc-1; ++i) {
    if (util::eq(argv[i], "-v", "--verbose")) verbose = true;
    if (util::eq(argv[i], "-t", "--tol")) {
      expect_another_arg(i, argc);
      ++i;
      tol = std::atof(argv[i]);
    }
  }

  // Always decorate file names with precision info
  std::string file1_fn(argv[argc-2]), file2_fn(argv[argc-1]);
  file1_fn += std::to_string(sizeof(Real));
  file2_fn += std::to_string(sizeof(Real));

  int rv = compare_files(file1_fn, file2_fn, tol, verbose);

  if (rv == 0) {
    std::cout << file1_fn << " and " << file2_fn << " appear to match" << std::endl;
  }
  else {
    std::cout << file1_fn << " and " << file2_fn << " do not match" << std::endl;
  }

  return rv == 0 ? 0 : 1;
}
