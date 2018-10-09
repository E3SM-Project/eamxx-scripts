#include "micro_sed_packnoi_kokkos.hpp"
#include "types.hpp"
#include "util.hpp"
#include "micro_sed_wrap.hpp"

int main(int argc, char** argv)
{
  common_main("micro_sed_packnoi_kokkos");

  Kokkos::initialize(argc, argv); {
    p3::micro_sed::micro_sed_func_kokkos_wrap<Real, p3::micro_sed::MicroSedFuncPackNoiKokkos<Real> >(ni, nk, dt, ts, kdir, repeat);
  } Kokkos::finalize();

  return 0;
}
