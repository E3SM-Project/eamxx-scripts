#include "micro_sed_vanilla.hpp"

#include <stdexcept>
#include <string>

typedef double Real;

int main(int argc, char** argv)
{
  if (argc != 8) {
    throw std::runtime_error(std::string("Usage: micro_sed kts kte ni nk its ite dt"));
  }

  int kts(atoi(argv[1])),
      kte(atoi(argv[2])),
      ni(atoi(argv[3])),
      nk(atoi(argv[4])),
      its(atoi(argv[5])),
      ite(atoi(argv[6]));


  Real dt(atof(argv[7]));


  p3::micro_sed_vanilla::micro_sed_func_vanilla_wrap<Real>(kts, kte, ni, nk, its, ite, dt);

  return 0;
}
