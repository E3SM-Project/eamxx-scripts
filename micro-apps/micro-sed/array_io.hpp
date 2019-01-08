#ifndef INCLUDE_ARRAY_IO_HPP
#define INCLUDE_ARRAY_IO_HPP

#include "types.hpp"

/*
 * This file contains interfaces for some simple routines for reading and
 * writing data arrays to files. They are compatible with Fortran and are used
 * in the micro-app, e.g., in p3_init for storing tables, and in the micro-app
 * drivers for recording rain-sed results for subsequent comparison.
 */

extern "C" {

// The following 3 routines work together, e.g., to optionally read in table
// data if an appropriate file exists, and otherwise to write table data after
// computing it for future use.

// Does filename exist as a file?
bool array_io_file_exists(const char* filename);
// Write the array a, having n elements, to file filename.
bool array_io_write(const char* filename, Real** a, const int n);
// Read from file filename to array a, expecting n elements.
bool array_io_read(const char* filename, Real** a, const int n);

// Micro-app driver routine to dump rain-sed data.
bool dump_all(const char* filename,
              const Real** qr, const Real** nr, const Real** th, const Real** dzq, const Real** pres, const Real** prt_liq,
              const int ni, const int nk, const Real dt, const int ts);

}

#endif
