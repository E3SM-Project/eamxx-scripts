#include "types.hpp"
#include "li_common.hpp"

extern "C" {

// This will link to the fortran reference implementation
void linear_interp_wrap(const int* km1, const int* km2, const int* ncol, const Real* minthresh, const int* repeat);

};

/*
 * This is the exe driver for the fortran reference rain-sed implementation.
 */

int main(int argc, char** argv)
{
  common_main("li_ref");

  linear_interp_wrap(&km1, &km2, &ncol, &minthresh, &repeat);

  return 0;
}
