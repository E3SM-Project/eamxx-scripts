# micro-sed
--------------------------------------------------
This directory contains the code needed to run/test the micro-sed micro-app
in standalone mode. The actual micro app code which will eventually be migrated
to p3 is in the ../p3/ directory.

## Implementations:
We created three implementations of rain sedimentation in this micro-app:

"ref" : The reference fortran implementation taken from p3. We took as little
        from p3 as we could and also tried to minimize modifications. The main
        modification was to add openmp threading so that we could get meaningful
        comparisons with C++ Kokkos implementations.
"vanilla" : This implementation was made by porting "ref" to C++ using Kokkos Views
            as the MD-array data structure. We made no effort to do any improvements
            besides the language change, so this is a direct, naive translation of the
            fortran.
"final" : This implementation represents what we think a good C++ rain sedimentaion
          should look like. This implementation should serve as a guideline/example for
          future C++/Kokkos rewrites of scientific fortran code. NOTE: most "final"
	  files are in ../p3/ instead of this directory because they will be migrated
	  to p3.

NOTE: renaming, adding, or removing implementations will have an impact outside this
directory, specifically on the testing and performance gathering convenience scripts.

For each implementation above (not including "ref" which just dives into fortran), you will see
(in either the micro-sed or p3 directory):
* A driver cpp file called p3_${impl}_driver.cpp. This gets compiled into an exe and actually
  executes the micro app. As such, it does the basic setup and argument parsing.
* A p3_${impl}.hpp file which has the declaration of the MicroSedFunc${impl} class which
  encapsulates that particular impl.
* A p3_${impl}_impl.hpp file, this has the definitions for the above class
* A p3_${impl}.cpp which has the explicit template instatiation for the above class.

## How the driver works:
--------------------------------------------------
* ../test_all calls the exe for each impl in order to do a timing comparison.
* That exe comes from p3_${impl}_driver.cpp, which calls common_main (which initializes
  stuff like dt and tables) and micro_sed_func_kokkos_wrap (which is templated on impl type).
  Both common_main and micro_sed_func_kokkos_wrap are from micro_app_common.hpp.
* micro_sed_func_kokkos_wrap initializes T, rain profile, etc, then calls micro_sed func inside
  a loop over the number of times the user wants to repeat the calc (for more robust timing)
  and handles timing and output stuff.
* micro_app_common.hpp accesses micro_sed_func by including p3_final.hpp

Discussion of micro_sed_func differs a bit between implementations and is described
in ../p3/README.md.


  

