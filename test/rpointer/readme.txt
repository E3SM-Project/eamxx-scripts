This directory contains a script, test-rpointer-mgr.sh, with inputs input*sh, to
test the rpointer manager against artificial crashes in the restart-write phase.

To run a test, edit an input file to point to your repo, 'git apply
rpointer_test.patch' to the repo (assuming PR 6248 is in your repo), and then run
    bash test-rpointer-mgr.sh input...sh all
This will do a sequence of runs:
    1. Run, write restarts, then crash in the second restart write.
    2. Continue twice.
    3. Run a monolithic case.
The script finishes when two logs are available with the final time stamp and
BFB information.

This procedure works for both the E3SM and SCREAM repos. Input files with
'scream' in the name must be run with the SCREAM repo.
