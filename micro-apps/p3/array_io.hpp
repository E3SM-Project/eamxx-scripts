#ifndef INCLUDE_ARRAY_IO_HPP
#define INCLUDE_ARRAY_IO_HPP

#include "types.hpp"

/*
 * This file contains interfaces for some simple routines for reading
 * and writing data arrays to files. They are intended to be
 * compatible with fortran and are used in the micro-app for storing
 * tables and recording rain-sed results from the fortran baseline
 * code.
 */

extern "C" {

bool array_io_file_exists(const char* filename);
bool array_io_write(const char* filename, Real** a, const int n);
bool array_io_read(const char* filename, Real** a, const int n);
bool dump_all(const char* filename,
              const Real** qr, const Real** nr, const Real** th, const Real** dzq, const Real** pres, const Real** prt_liq,
              const int ni, const int nk, const Real dt, const int ts);

}

#endif
