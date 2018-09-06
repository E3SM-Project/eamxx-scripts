#include "types.hpp"
#include "util.hpp"
#include "initial_conditions.hpp"

#include <vector>
#include <iostream>
#include <exception>

extern "C" {

// This is for Fortran only, so pack in Fortran order.
void fully_populate_input_data(Int ni, Int nk, Int kdir, Real** qr, Real** nr, Real** th, Real** dzq, Real** pres)
{
  ic::MicroSedData<Real> data(ni, nk);
  populate(data, kdir);

  for (auto item : { std::make_pair(data.qr, *qr), std::make_pair(data.nr, *nr), std::make_pair(data.th, *th),
                     std::make_pair(data.dzq, *dzq), std::make_pair(data.pres, *pres) }) {
    for (int k = 0; k < nk; ++k)
      for (int i = 0; i < ni; ++i)
        item.second[ni*k + i] = item.first[k];
  }
}

}
