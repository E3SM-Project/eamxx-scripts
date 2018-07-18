#ifndef INCLUDE_ARRAY_IO_HPP
#define INCLUDE_ARRAY_IO_HPP

namespace array_io {
typedef float freal;
}

extern "C" {
  bool array_io_file_exists(const char* filename);
  bool array_io_write(const char* filename, array_io::freal** a, const int n);
  bool array_io_read(const char* filename, array_io::freal** a, const int n);
}

#endif
