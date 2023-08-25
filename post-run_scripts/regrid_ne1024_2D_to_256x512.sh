#!/bin/bash                                 
#$ -cwd                                     

# regrids ne1024np4 output to 70 km fv grid 
# Tested script having requested 4 interactive nodes on cori-haswell                                          
# Mainly for regridding 2D variables        

# ============================================================                                                

# Location of the model output              
# USER needs to change this                 
fileloc='/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/20201112.SCREAMv0dyamond2.F2010-SCREAM-HR-DYAMOND2.ne1024pg2_r0125_oRRS18to6v3.cori-knl.1536x8x16/run/'

# Location of the regridded files           
# USER needs to change this
outputloc='/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201112/regridded/'

# =============================================================                                               

# Calls a bash script that retrieves, regrids, and concatenates each variable data                            
# Variables can be added or removed based on the interest of the USER, but currently regrids h0, h1, h2, h3, h4 files in order                                                  

# USER: Compare the fincl list with the list below to make sure they correspond                               
h0_files=("CLDLOW" "CLDMED" "CLDHGH" "CLDTOT" "TMCLDLIQ" "TMCLDICE" "TMRAINQM" "TMCLDRIM" "TMQ" "CAPE" "CIN")
h1_files=("PS" "TS" "TREFHT" "QREFHT" "PRECT" "PRESL" "WINDSPD_10M" "TAUX" "TAUY" "SHFLX" "LHFLX")
h2_files=("FSNTOA" "FLNT" "FLNTC" "FSNTOAC" "FSNS" "FSDS" "FLNS" "FLDS")
h3_files=("RH200" "RH500" "RH700" "RH850" "OMEGA200" "OMEGA500" "OMEGA700" "OMEGA850" "Z200" "Z500" "Z700" "Z850")
h4_files=("PS" "PSL" "TMNUMLIQ" "TMNUMICE" "TMNUMRAI")
echo "  "
echo " ---------------- Starting h0 files -------------------- "
echo "  "

for var in ${h0_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h0 $outputloc
done

echo "  "
echo " ---------------- Starting h1 files -------------------- "
echo "  "

for var in ${h1_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h1 $outputloc
done

echo "  "
echo " ---------------- Starting h2 files -------------------- "
echo "  "


for var in ${h2_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h2 $outputloc
done

echo "  "
echo " ---------------- Starting h3 files -------------------- "
echo "  "


for var in ${h3_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h3 $outputloc
done

echo "  "
echo " ---------------- Starting h4 files -------------------- "
echo "  "

for var in ${h4_files[*]}; do
    echo " "
    echo $var
    bash regrid_multi_ne1024_files_lowrezOutput.sh $fileloc $var h4 $outputloc
done


echo "  "
echo "  "
echo "  "
echo "   DONE!  "
echo "  "
echo "  "
