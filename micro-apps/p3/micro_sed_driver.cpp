#include "micro_sed_vanilla.hpp"
#include "types.hpp"
#include "util.hpp"

extern "C" {

void micro_sed_func_wrap(const int* ni, const int* nk, const Real* dt, const int* ts, const int* kdir);

};

int main(int argc, char** argv)
{
  common_main("micro_sed");

  micro_sed_func_wrap(&ni, &nk, &dt, &ts, &kdir);

  return 0;
}
