#include "micro_sed_vanilla_kokkos.hpp"
#include "types.hpp"
#include "util.hpp"

int main(int argc, char** argv)
{
  micro_throw_if(argc != 6, "Usage: micro_sed_vanilla_kokkos ni nk time_step_len num_steps kdir");

  int ni(atoi(argv[1])),
      nk(atoi(argv[2])),
      ts(atoi(argv[4])),
      kdir(atoi(argv[5]));

  Real dt(atof(argv[3]));

  micro_throw_if(kdir != -1 && kdir != 1, "kdir must be -1 or 1");

  Kokkos::initialize(argc, argv); {
    p3::micro_sed_vanilla::p3_init_cpp<Real>();
    p3::micro_sed_vanilla::micro_sed_func_vanilla_kokkos_wrap<Real>(ni, nk, dt, ts, kdir);
  } Kokkos::finalize();

  return 0;
}
