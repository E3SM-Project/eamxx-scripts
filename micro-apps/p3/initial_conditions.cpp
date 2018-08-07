#include "types.hpp"
#include "util.hpp"
#include "initial_conditions.hpp"

#include <vector>
#include <iostream>
#include <exception>

extern "C" {

void fully_populate_input_data(Int ni, Int nk, Real** qr, Real** nr, Real** th, Real** dzq, Real** pres)
{
  ic::MicroSedData<Real> data(ni, nk);
  populate(data);

  for (auto item : { std::make_pair(data.qr, *qr), std::make_pair(data.nr, *nr), std::make_pair(data.th, *th),
                     std::make_pair(data.dzq, *dzq), std::make_pair(data.pres, *pres) }) {
    std::copy(item.first, item.first + ni*nk, item.second);
  }
}

}
