#include "array_io.hpp"
#include "util.hpp"

#include <sys/stat.h>

#include <iostream>

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
}
