# p3
---------------------------
This directory contains p3 code from micro-apps which we intend to migrate to
(the scream fork of) the E3SM repo. In its current location, this code is used
by driver code which is in ../micro-sed/.

## Code Structure:
--------------------------
* p3_final.hpp declares/points to (rather than defines) the MicroSedFinalKokkos class
  which contains micro_sed_func. If using GPUs (which require template definitions to 
  be available for each translation unit), Micro_SedFinalKokkos also defines micro_sed_func
  by including pe_final_impl.hpp. Otherwise ETI is used. 
* p3_final.cpp instantiates MicroSedKokkos and defines it by including p3_final_impl.hpp. 
* p3_final_impl.hpp is the main code implementing the micro app. Note that it uses the Functions
  class included via p3_functions.hpp which itself uses find, upwind, and table functions
  handled via ETI.



