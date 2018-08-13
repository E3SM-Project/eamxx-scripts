#include "micro_sed_vanilla_kokkos.hpp"
#include "types.hpp"
#include "util.hpp"

int main(int argc, char** argv)
{
  micro_throw_if(argc != 9, "Usage: micro_sed kts kte ni nk its ite dt");

  int kts(atoi(argv[1])),
      kte(atoi(argv[2])),
      ni(atoi(argv[3])),
      nk(atoi(argv[4])),
      its(atoi(argv[5])),
      ite(atoi(argv[6])),
      ts(atoi(argv[8]));

  Real dt(atof(argv[7]));

  p3::micro_sed_vanilla::p3_init_cpp_kokkos<Real>();

  Kokkos::initialize(argc, argv); {
    p3::micro_sed_vanilla::micro_sed_func_vanilla_kokkos_wrap<Real>(kts, kte, ni, nk, its, ite, dt, ts);
  } Kokkos::finalize_all();

  return 0;
}
