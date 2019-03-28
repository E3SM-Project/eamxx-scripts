#include "li_common.hpp"
#include "../micro-sed/cmp.hpp"

extern "C" {

void populate_li_input_from_fortran(int km1, int km2, Real** x1_i, Real** y1_i, Real** x2_i)
{
  li::populate_li_input(km1, km2, *x1_i, *y1_i, *x2_i);
}

bool dump_all_li(const char* filename,
                 const Real** y2,
                 const int ncol, const int km1, const int km2, const Real minthresh)
{
  const int size = ncol * km2;
  std::vector<Real> y2_cpp(size);
  cmp::transpose<cmp::TransposeDirection::f2c>(*y2, y2_cpp.data(), ncol, km2);
  li::dump_to_file_li(filename, y2_cpp.data(), ncol, km1, km2, minthresh);
  return true;
}

}

namespace li {



} // namespace li
