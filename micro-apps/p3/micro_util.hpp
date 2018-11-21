#ifndef MICRO_UTIL_HPP
#define MICRO_UTIL_HPP

#include "scream_assert.hpp"
#include "scream_arch.hpp"

/*
 * Utilities specific to the micro app.
 *
 * This stuff probably won't be needed in the main scream repo.
 */

#define common_main(exename)                                            \
  util::initialize();                                                   \
  micro_require_msg(argc == 7, "Usage: " #exename " ni nk time_step_len num_steps kdir repeat"); \
  int ni(atoi(argv[1])), nk(atoi(argv[2])), ts(atoi(argv[4])), kdir(atoi(argv[5])), repeat(atoi(argv[6])); \
  Real dt(atof(argv[3]));                                               \
  micro_require_msg(kdir == -1 || kdir == 1, "kdir must be -1 or 1");   \
  p3::micro_sed::p3_init_cpp<Real>()

#endif
