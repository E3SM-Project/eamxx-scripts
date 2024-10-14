# eamxx-scripts
--------------------------------------------------
NOTE: many scripts in this repo were created for previous scream configurations and simulations.
And none of the materials in this repo are curated or supported.

This is where scripts related to the Simple Cloud Resolving
E3SM Atmosphere Model live. This is also where micro-app code and related testing
scripts are put before migrating finished code to the (scream fork of the) E3SM repo.


## Top-level Code structure (as of 10/14/24):
--------------------------------------------------
* f2py/p3/:		 Initial standalone version of p3 which Peter created using f2py.
* perf-scripts/: 	         Stuff for making performance plots
* post-run_scripts/: 	     Includes some scripts used for checking and post-processing scream output.
* preprocessing_scripts/:    Includes scripts for creating initial conditions for eamxx
* micro-apps/:   	         Includes code for each micro-app we've generated
* run_scripts/:              Includes run scripts previously used to run scream production simulations. 
* util/:                     Includes job monitor script and script to generate a map for regional outputs.
* v1_output/:                Includes sample output .yaml files used in previous simulations.