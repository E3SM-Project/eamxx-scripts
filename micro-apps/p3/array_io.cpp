#include "array_io.hpp"
#include "util.hpp"
#include "micro_sed_vanilla.hpp"
#include "cmp.hpp"

#include <sys/stat.h>

#include <iostream>
#include <vector>

namespace array_io {

template <typename Scalar>
void write (const char* filename, Scalar* a, const int n) {
  util::FILEPtr fid(fopen(filename, "w"));
  micro_throw_if( ! fid, "Could not open " << filename << " for writing.");
  util::write<int>(&n, 1, fid);
  util::write<Scalar>(a, n, fid);
}

template <typename Scalar>
void read (const char* filename, Scalar* a, const int n) {
  util::FILEPtr fid(fopen(filename, "r"));
  micro_throw_if( ! fid, "Could not open " << filename << " for reading.");
  int n_file;
  util::read<int>(&n_file, 1, fid);
  micro_throw_if(n_file != n, "Expected " << n << " but got " << n_file);
  util::read<Scalar>(a, n, fid);
}

} // namespace array_io

extern "C" {

bool array_io_file_exists (const char* filename) {
  struct stat s;
  const bool exists = stat(filename, &s) == 0;
  return exists;
}

bool array_io_write (const char* filename, Real** a, const int n) {
  try {
    array_io::write(filename, *a, n);
    return true;
  } catch (std::exception& e) {
    std::cerr << "array_io_write failed with: " << e.what() << "\n";
    return false;
  }
}

bool array_io_read (const char* filename, Real** a, const int n) {
  try {
    array_io::read(filename, *a, n);
    return true;
  } catch (std::exception& e) {
    std::cerr << "array_io_read failed with: " << e.what() << "\n";
    return false;
  }
}

bool dump_all(const char* filename,
              const Real** qr, const Real** nr, const Real** th, const Real** dzq, const Real** pres, const Real** prt_liq,
              const int ni, const int nk, const Real dt, const int ts)
{
  try {
    const int size = ni * nk;
    std::vector<Real> qr_t(size), nr_t(size), th_t(size), dzq_t(size), pres_t(size);
    for (auto item : { std::make_pair(*qr, qr_t.data()), std::make_pair(*nr, nr_t.data()), std::make_pair(*th, th_t.data()),
          std::make_pair(*dzq, dzq_t.data()), std::make_pair(*pres, pres_t.data()) }) {
      cmp::transpose<cmp::TransposeDirection::f2c>(item.first, item.second, ni, nk);
    }

    p3::micro_sed_vanilla::dump_to_file(filename, qr_t.data(), nr_t.data(), th_t.data(), dzq_t.data(), pres_t.data(),
                                        *prt_liq, ni, nk, dt, ts);
    return true;
  } catch (std::exception& e) {
    std::cerr << "dump_all failed with: " << e.what() << "\n";
    return false;
  }
}


}
