#include "li_common.hpp"

extern "C" {

void populate_li_input_from_fortran(int km1, int km2, Real** x1_i, Real** y1_i, Real** x2_i)
{
  li::populate_li_input(km1, km2, *x1_i, *y1_i, *x2_i);
}

bool dump_all_li(const char* filename,
                 const Real** y2,
                 const int ncol, const int km1, const int km2, const Real minthresh)
{
  return true;
}

}

namespace li {



} // namespace li
