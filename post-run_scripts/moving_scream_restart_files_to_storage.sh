#!/bin/bash
#$ -cwd

# Script moves data to HPSS and verifies the integrity of the files

# ============================================================

# Location of the model output
# USER needs to change this
fileloc='/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SCREAMv0.SCREAM-DY2.ne1024pg2.Feb11Case.20210115/run/'

fileprefix='SCREAMv0.SCREAM-DY2.ne1024pg2.Feb11Case.20210115'
# Location of the regridded files
# USER needs to change this
hpssloc='/home/t/terai/Production_runs/SCREAMv0.SCREAM-DY2.ne1024pg2.Feb11Case.20210115/'

# =============================================================


echo "  "
echo " ---------------- Movement of data  -------------------- "
echo "  "

cd ${fileloc}

for f in ${fileprefix}.*.r*2020-02-13-00000* ; do
    
    echo "moving $f"
    hsi -P "cd ${hpssloc}; cput -c on $f" >> moving_restart_log_Feb11.txt
done

#for f in ${fileprefix}.c*.r*.2020-03-01-00000* ; do
#    
#    echo "moving $f"
#    hsi -P "cd ${hpssloc}; cput -c on $f" >> moving_restart_log.txt
#done

for f in ${fileprefix}.*.r*2020-02-13-00000* ; do

    echo "verifying $f"
    hsi -P hashverify ${hpssloc}${f} >> moving_verify_restart_Feb11.txt
done

#for f in ${fileprefix}.c*.r*.2020-03-01-00000* ; do
#
#    echo "verifying $f"
#    hsi -P hashverify ${hpssloc}${f} >> moving_verify_restart.txt
#done


echo "  "
echo "  "

echo "  "
echo "  "
echo "  "
echo "   DONE!  "
echo "  "
echo "  "
echo "  ******** "  
