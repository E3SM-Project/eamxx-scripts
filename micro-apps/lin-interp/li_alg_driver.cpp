#include "li_alg.hpp"
#include "types.hpp"
#include "li_common.hpp"

/*
 * This is the exe driver for the vanilla kokkos rain-sed implementation.
 */

int main(int argc, char** argv)
{
  common_main("li_alg");

  Kokkos::initialize(argc, argv); {
    li::lin_interp_func_wrap<Real, li::LiAlg<Real> >(ncol, km1, km2, minthresh, repeat);
  } Kokkos::finalize();

  return 0;
}
