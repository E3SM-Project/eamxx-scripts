#!/bin/bash
#$ -cwd

# Script moves data to HPSS and verifies the integrity of the files

# ============================================================

# Location of the model output
# USER needs to change this
fileloc='/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/20201112.SCREAMv0dyamond2.F2010-SCREAM-HR-DYAMOND2.ne1024pg2_r0125_oRRS18to6v3.cori-knl.1536x8x16/run/'

fileprefix='20201112.SCREAMv0dyamond2.F2010-SCREAM-HR-DYAMOND2.ne1024pg2_r0125_oRRS18to6v3.cori-knl.1536x8x16'
# Location of the regridded files
# USER needs to change this
hpssloc='/home/t/terai/Production_runs/20201112.SCREAMv0dyamond2.F2010-SCREAM-HR-DYAMOND2.ne1024pg2_r0125_oRRS18to6v3.cori-knl.1536x8x16/'

# =============================================================


echo "  "
echo " ---------------- Movement of data  -------------------- "
echo "  "

cd ${fileloc}

for f in ${fileprefix}*.r*.* ; do
    
    echo "moving $f"
    time hsi "cd ${hpssloc}; put -c on $f"
done

for f in ${fileprefix}*.r*.* ; do

    echo "verifying $f"
    time hsi hashverify ${hpssloc}${f}
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
