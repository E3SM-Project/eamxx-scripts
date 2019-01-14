# p3
---------------------------
This directory contains p3 code from micro-apps which we intend to migrate to
(the scream fork of) the E3SM repo. In its current location, this code is used
by driver code which is in ../micro-sed/.

## Code Structure:
--------------------------

* p3_final.hpp declares/points to (rather than defines) the MicroSedFinalKokkos
  class which contains micro_sed_func. If using GPUs, Micro_SedFinalKokkos also
  defines micro_sed_func by including p3_final_impl.hpp, thus making all
  definitions available to the kernel. Otherwise, explicit template
  instantiation (ETI) is used.
* In general, ETI separates the build into multiple small translation units to
  improve build time. Each translation unit can define template classes and
  functions for specific template parameters, and these can then be used by
  other translation units at link time. ETI (and, more generally, placing
  definitions, whether templated or not, in many small translation units) is
  helpful on CPU/KNL so that the Intel compiler can spend time in each small
  translation unit auto-vectorizing, as guided by the Pack mechanic. The GPU
  requires either that function definitions be available to each kernel, or that
  the relocatable device code (RDC) flag be used. In the first case, ETI is thus
  not permitted. The second case may be an option in the future, but it has
  performance downsides. Since nvcc does not auto-vectorize the code (GPU uses a
  different mechanic), the build time is currently quite tolerable without RDC.
* p3_final.cpp instantiates MicroSedKokkos and defines it by including p3_final_impl.hpp. 
* p3_final_impl.hpp is the main code implementing the micro app. Note that it uses the Functions
  class included via p3_functions.hpp which itself uses find, upwind, and table functions
  handled via ETI.



