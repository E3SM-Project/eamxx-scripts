#include "p3_common.hpp"
#include "types.hpp"
#include "micro_util.hpp"

extern "C" {

void micro_sed_func_wrap(const int* ni, const int* nk, const Real* dt, const int* ts, const int* kdir, const int* repeat);

};

int main(int argc, char** argv)
{
  common_main("p3_ref");

  micro_sed_func_wrap(&ni, &nk, &dt, &ts, &kdir, &repeat);

  return 0;
}
