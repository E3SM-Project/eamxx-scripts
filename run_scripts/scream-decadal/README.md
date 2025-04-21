Runscripts used for the SCREAM decadal AMIP-style run. There are three runscripts used:

run.decadal-amip.sh: used for setting up the original run (the main runscript)

run.decadal-amip-alternate.sh: used for setting up a pseudo-branch-run with the alternative fortran compiler to evaluate whether this configuration would be more stable. This parallel effort was eventually abandoned, as issues with the original run seemed to resolve.

run.decadal-amip-fill-missing.sh: used for setting up pseudo-branch-runs to repeat sections of the run with missing or corrupt data arising due to run crashes
