#!/bin/bash
#$ -cwd

# regrids ne1024np4 output to 70 km fv grid
# Tested script having requested 4 interactive nodes on cori-haswell 
# Mainly for regridding 2D variables

# ============================================================

# Location of the model output
# USER needs to change this
fileloc='/global/cscratch1/sd/terai/E3SM_simulations/master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1536x8x16.DY2_Nov05branch_SHOC_P3_AB_bugfix.20201105-16/run/'

# Location of the regridded files
# user needs to change this
outputloc='/global/cscratch1/sd/terai/SCREAM/output/Prelim_runs/master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.20201105-16/'

# =============================================================

# Calls a bash script that retrieves, regrids, and concatenates each variable data
# Variables can be added or removed based on the interest of the USER, but currently regrids h0, h1, h2, h3, h4 files in order

h5_files=("PS")
h6_files=("U" "V" "OMEGA")
h7_files=("T" "Q")
h8_files=("CLDLIQ" "CLDICE")
h9_files=("TOT_ICLD_VISTAU" "EMIS" "CLOUD")

echo "  "
echo " ---------------- Starting h5 files -------------------- "
echo "  "

for var in ${h5_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h5 $outputloc
done

echo "  "
echo " ---------------- Starting h6 files -------------------- "
echo "  "

for var in ${h6_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h6 $outputloc
done

echo "  "
echo " ---------------- Starting h7 files -------------------- "
echo "  "


for var in ${h7_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h7 $outputloc
done

echo "  "
echo " ---------------- Starting h8 files -------------------- "
echo "  "


for var in ${h8_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h8 $outputloc
done

echo "  "
echo " ---------------- Starting h9 files -------------------- "
echo "  "

for var in ${h9_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h9 $outputloc
done



echo "  "
echo "  "
echo "  "
echo "   DONE!  "
echo "  "
echo "  "
echo "  ******** "  
