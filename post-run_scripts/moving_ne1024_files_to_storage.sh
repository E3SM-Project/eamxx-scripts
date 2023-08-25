#!/bin/bash
#$ -cwd

# Script moves data to HPSS and verifies the integrity of the files

# ============================================================

# Location of the model output
# USER needs to change this
fileloc='/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SCREAMv0.SCREAM-DY2.ne1024pg2.Feb11Case.20210115/run/'

# Location of the regridded files
# USER needs to change this
hpssloc='/home/t/terai/Production_runs/SCREAMv0.SCREAM-DY2.ne1024pg2.Feb11Case.20210115/'

# =============================================================


echo "  "
echo " ---------------- Movement of data  -------------------- "
echo "  "

cd ${fileloc}

for f in SCREAMv0.SCREAM-DY2.ne1024pg2.Feb11Case.20210115.eam.h*.nc ; do
    
    echo "moving $f"
    hsi "cd ${hpssloc}; put -c on $f"
done

for f in SCREAMv0.SCREAM-DY2.ne1024pg2.Feb11Case.20210115.eam.h*.nc ; do

    echo "verifying $f"
    hsi -P hashverify ${hpssloc}${f} >> move_verify_Feb11.txt
done



echo "  "
echo "  "

echo "  "
echo "  "
echo "  "
echo "   DONE!  "
echo "  "
echo "  "
echo "  ******** "  
