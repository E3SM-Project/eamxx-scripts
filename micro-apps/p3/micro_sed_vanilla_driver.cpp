#include "micro_sed_vanilla.hpp"
#include "util.hpp"

#include <stdexcept>
#include <string>

using Real = double;

int main(int argc, char** argv)
{
  micro_throw_if(argc != 8, "Usage: micro_sed kts kte ni nk its ite dt");

  int kts(atoi(argv[1])),
      kte(atoi(argv[2])),
      ni(atoi(argv[3])),
      nk(atoi(argv[4])),
      its(atoi(argv[5])),
      ite(atoi(argv[6]));

  Real dt(atof(argv[7]));

  p3::micro_sed_vanilla::p3_init_cpp<Real>();

  p3::micro_sed_vanilla::micro_sed_func_vanilla_wrap<Real>(kts, kte, ni, nk, its, ite, dt);

  return 0;
}
