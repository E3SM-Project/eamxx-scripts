#include "micro_sed_workspace_kokkos.hpp"
#include "types.hpp"
#include "util.hpp"

int main(int argc, char** argv)
{
  common_main("micro_sed_workspace_kokkos");

  Kokkos::initialize(argc, argv); {
    p3::micro_sed_vanilla::micro_sed_func_workspace_kokkos_wrap<Real>(ni, nk, dt, ts, kdir);
  } Kokkos::finalize();

  return 0;
}
