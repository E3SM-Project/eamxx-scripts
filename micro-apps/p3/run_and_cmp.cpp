#include "types.hpp"
#include "util.hpp"

extern "C" void
micro_sed_func_c(
  Int kts, Int kte, Int ni, Int nk, Int its, Int ite, const Real* dt,
  Real* qr, Real* nr, const Real* th, const Real* dzq, const Real* pres,
  Real* prt_liq);

int main (int argc, char** argv) {
  
}
