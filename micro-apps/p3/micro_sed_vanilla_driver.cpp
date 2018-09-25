#include "micro_sed_vanilla.hpp"
#include "types.hpp"
#include "util.hpp"
#include "micro_sed_wrap.hpp"

int main(int argc, char** argv)
{
  common_main("micro_sed_vanilla");

  p3::micro_sed::micro_sed_func_vanilla_wrap<Real>(ni, nk, dt, ts, kdir, repeat);

  return 0;
}
