#include "types.hpp"
#include "util.hpp"
#include "initial_conditions.hpp"

#include <vector>
#include <iostream>
#include <exception>

extern "C" {

// This is for Fortran only, so pack in Fortran order.
void populate_input_from_fortran(Int nk, Int kdir, Real** qr, Real** nr, Real** th, Real** dzq, Real** pres)
{
  ic::MicroSedData<Real> data(1, nk);
  populate(data, kdir);

  for (auto item : { std::make_pair(data.qr, *qr), std::make_pair(data.nr, *nr), std::make_pair(data.th, *th),
                     std::make_pair(data.dzq, *dzq), std::make_pair(data.pres, *pres) }) {
    for (int k = 0; k < nk; ++k)
      item.second[k] = item.first[k];
  }
}

}
