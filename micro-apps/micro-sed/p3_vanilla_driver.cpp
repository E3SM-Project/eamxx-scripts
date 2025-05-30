#include "p3_vanilla.hpp"
#include "types.hpp"
#include "micro_app_common.hpp"

/*
 * This is the exe driver for the vanilla kokkos rain-sed implementation.
 */

int main(int argc, char** argv)
{
  common_main("p3_vanilla");

  Kokkos::initialize(argc, argv); {
    p3::micro_sed::micro_sed_func_kokkos_wrap<Real, p3::micro_sed::MicroSedFuncVanillaKokkos<Real> >(ni, nk, dt, ts, kdir, repeat);
  } Kokkos::finalize();

  return 0;
}
