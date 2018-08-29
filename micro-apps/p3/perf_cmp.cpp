#include "types.hpp"
#include "util.hpp"
#include "cmp.hpp"

#include <vector>
#include <iostream>
#include <exception>

static void expect_another_arg (Int i, Int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}

template <typename Scalar>
Scalar compare_property(const std::string& file1_fn, const std::string& file2_fn,
                        const util::FILEPtr& fid1, const util::FILEPtr& fid2, const char* name,
                        bool verbose)
{
  Scalar item1, item2;

  util::read(&item1, 1, fid1);
  util::read(&item2, 1, fid2);

  micro_throw_if(item1 != item2, "Files " << file1_fn << " and " << file2_fn <<
                 " have mismatch in basic property '" << name << "', " << item1 << " != " << item2);

  if (verbose) std::cout << name << " matches, = " << item1 << std::endl;

  return item1;
}

int compare_files(const std::string& file1_fn, const std::string& file2_fn, Real tol, bool verbose)
{

  util::FILEPtr fid1(fopen(file1_fn.c_str(), "r")), fid2(fopen(file2_fn.c_str(), "r"));

  int ni = compare_property<int>(file1_fn, file2_fn, fid1, fid2, "ni", verbose);
  int nk = compare_property<int>(file1_fn, file2_fn, fid1, fid2, "nk", verbose);
  compare_property<Real>(file1_fn, file2_fn, fid1, fid2, "dt", verbose);
  compare_property<int>(file1_fn, file2_fn, fid1, fid2, "ts", verbose);

  const int size = ni * nk;
  int bad_cnt = 0;
  const char* prt_liq = "prt_liq";
  for (auto item : { "qr", "nr", "th", "dzq", "pres", prt_liq } ) {
    const int curr_size = item == prt_liq ? ni : size;
    std::vector<Real> vec1(curr_size), vec2(curr_size);
    util::read(vec1.data(), curr_size, fid1);
    util::read(vec2.data(), curr_size, fid2);
    int curr_bad = cmp::compare(item, vec1.data(), vec2.data(), curr_size, tol, verbose);
    bad_cnt += curr_bad;
    if (curr_bad == 0 && verbose) std::cout << item << " matched " << curr_size << " items." << std::endl;
  }

  return bad_cnt;
}

int main (int argc, char** argv) {
  if (argc < 3) {
    std::cout <<
      argv[0] << " [options] file1 file2\n"
      "Options:\n"
      "  -v        Run with extra verbose output.\n"
      "  -t <tol>  Tolerance for relative error.\n";
    return -1;
  }

  Real tol = 0;
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

  if (util::is_single_precision<Real>::value &&
      (file1_fn.find("fortran") != file2_fn.find("fortran")) &&
      tol < util::TOL) {
    std::cout << "Due to single-precision and comparing fortran to non-fortran, bumping tol to " << util::TOL << std::endl;
    tol = util::TOL;
  }

  if ( (file1_fn.find("kokkos") != std::string::npos || file2_fn.find("kokkos") != std::string::npos) &&
       tol < util::TOL) {
    std::cout << "Due to presence of kokkos, bumping tol to " << util::TOL << std::endl;
    tol = util::TOL;
  }

  int rv = compare_files(file1_fn, file2_fn, tol, verbose);

  if (rv == 0) {
    std::cout << file1_fn << " and " << file2_fn << " appear to match" << std::endl;
  }
  else {
    std::cout << file1_fn << " and " << file2_fn << " do not match" << std::endl;
  }

  return rv;
}
