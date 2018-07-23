#include "types.hpp"
#include "util.hpp"

#include <iostream>

extern "C" void
micro_sed_func_c(
  Int kts, Int kte, Int ni, Int nk, Int its, Int ite, const Real* dt,
  Real* qr, Real* nr, const Real* th, const Real* dzq, const Real* pres,
  Real* prt_liq);

int main (int argc, char** argv) {
  if (argc == 0 || argc > 3) {
    std::cout << argv[0]
              << " -g baseline-filename-out: Generate baseline file.\n"
              << argv[0]
              << " baseline-filename-in: Run tests and compare against this"
              << " baseline file.\n";
  }
}
