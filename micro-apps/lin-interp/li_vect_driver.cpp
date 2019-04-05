#include "li_vect.hpp"
#include "types.hpp"
#include "li_common.hpp"

/*
 * This is the exe driver for the vanilla kokkos rain-sed implementation.
 */

int main(int argc, char** argv)
{
  common_main("li_vect");

  Kokkos::initialize(argc, argv); {
    li::lin_interp_func_wrap_kokkos<Real, li::LiVect<Real> >(ncol, km1, km2, minthresh, repeat);
  } Kokkos::finalize();

  return 0;
}
