#ifndef INCLUDE_ARRAY_IO_HPP
#define INCLUDE_ARRAY_IO_HPP

#include "types.hpp"

extern "C" {
  bool array_io_file_exists(const char* filename);
  bool array_io_write(const char* filename, Real** a, const int n);
  bool array_io_read(const char* filename, Real** a, const int n);
}

#endif
