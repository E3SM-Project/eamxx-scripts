# p3
==========================
This directory contains p3 code from micro-apps which we intend to migrate to
(the scream fork of) the E3SM repo. In its current location, this code is used
by driver code which is in ../micro-sed/.

## Code Structure:
--------------------------
* p3_final.hpp declares/points to (rather than defines) the MicroSedFinalKokkos class
  which contains micro_sed_func. If using GPUs, Micro_SedFinalKokkos also defines micro_sed_func
  by including pe_final_impl.hpp. Otherwise ETI is used. WHY IS ETI USED FOR NON-GPU ONLY?
* p3_final.cpp instantiates MicroSedKokkos and defines it by including p3_final_impl.hpp. I
  STILL DON'T UNDERSTAND THE CIRCUMSTANCES UNDER WHICH P3_FINAL.CPP WOULD BE CALLED.
* p3_final_impl.hpp is the main code implementing the micro app. Note that it uses the Functions
  class included via p3_functions.hpp which itself uses find, upwind, and table functions
  handled via ETI.



