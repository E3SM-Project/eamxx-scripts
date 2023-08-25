#!/bin/bash
#$ -cwd

# SCREAM history has variable long names that are cut off, script takes the long name from the first file and appends the correct name

# ============================================================

# Location of the model output
# USER needs to change this
#fileloc='/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127/'
fileloc='/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127/run/'

fileprefix='SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h9'

# =============================================================

# some input info

h0_files=("CLDLOW" "CLDMED" "CLDHGH" "CLDTOT" "TMCLDLIQ" "TMCLDICE" "TMRAINQM" "TMCLDRIM" "TMQ" "CAPE" "CIN")
h1_files=("PS" "TS" "TREFHT" "QREFHT" "PRECT" "PRECSL" "WINDSPD_10M" "TAUX" "TAUY" "SHFLX" "LHFLX")
h2_files=("FSNTOA" "FLNT" "FLNTC" "FSNTOAC" "FSNS" "FSDS" "FLNS" "FLDS")
h3_files=("RH200" "RH500" "RH700" "RH850" "OMEGA200" "OMEGA500" "OMEGA700" "OMEGA850" "Z200" "Z500" "Z700" "Z850")
h4_files=("PS" "PSL" "TMNUMLIQ" "TMNUMICE" "TMNUMRAI")
h5_files=("U" "V")
h6_files=("T" "Q")
h7_files=("CLDLIQ" "CLDICE")
h8_files=("CLOUD" "OMEGA")
h9_files=("EMIS" "TOT_ICLD_VISTAU")




# Following functions written 

ncvarlst() { ncks --trd -m ${1} | grep -E ': type' | cut -f 1 -d ' ' | sed 's/://' | sort ; }

ncattget() { ncks --trd -M -m ${3} | grep -E -i "^${2} attribute [0-9]+: ${1}" | cut -f 11- -d ' ' | sort ; }



echo " "
echo " --------------- Get variable longnames ------------------- "
echo " "

cd ${fileloc}

ncvarlst SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h9.2020-01-20-00000.nc

for varnamec in ${h9_files[*]} ; do
    echo " "
    echo $varnamec
    ncattget long_name $varnamec SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h9.2020-01-20-00000.nc
done


echo "  "
echo " ---------------- Changing attributes  -------------------- "
echo "  "

echo ${fileprefix}
for filenamec in ${fileprefix}*20-0[123]-*-00000.nc ; do
    echo ${filenamec}
    #ncatted -a long_name,CAPE,m,c,"Convectively available potential energy" -a long_name,CLDLOW,m,c,"Vertically-integrated low cloud" -along_name,CIN,m,c,"Convective inhibition" -a long_name,CLDHGH,m,c,"Vertically-integrated high cloud" -a long_name,CLDMED,m,c,"Vertically-integrated mid-level cloud" -a long_name,CLDTOT,m,c,"Vertically-integrated total cloud" -a long_name,TMCLDICE,m,c,"CLDICE column burden" -a long_name,TMCLDLIQ,m,c,"CLDLIQ column burden" -a long_name,TMCLDRIM,m,c,"CLDRIM column burden" -a long_name,TMQ,m,c,"Total (vertically integrated) precipitable water" -a long_name,TMRAINQM,m,c,"RAINQM column burden"  ${filenamec}
    #ncatted -a long_name,PS,m,c,"Surface pressure" -a long_name,TS,m,c,"Surface temperature (radiative)" -a long_name,TREFHT,m,c,"Reference height temperature" -a long_name,QREFHT,m,c,"Reference height humidity" -a long_name,PRECT,m,c,"Total (convective and large-scale) precipitation rate (liq + ice)" -a long_name,PRECSL,m,c,"Large-scale (stable) snow rate (water equivalent)" -a long_name,WINDSPD_10M,m,c,"10m wind speed" -a long_name,TAUX,m,c,"Zonal surface stress" -a long_name,TAUY,m,c,"Meridional surface stress" -a long_name,SHFLX,m,c,"Surface sensible heat flux" -a long_name,LHFLX,m,c,"Surface latent heat flux" ${filenamec}
    #ncatted -a long_name,FSNTOA,m,c,"Net solar flux at top of atmosphere" -a long_name,FLNT,m,c,"Net longwave flux at top of model" -a long_name,FLNTC,m,c,"Clearsky net longwave flux at top of model" -a long_name,FSNTOAC,m,c,"Clearsky net solar flux at top of atmosphere" -a long_name,FSNS,m,c,"Net solar flux at surface" -a long_name,FLNS,m,c,"Net longwave flux at surface" -a long_name,FSDS,m,c,"Downwelling solar flux at surface" -a long_name,FLDS,m,c,"Downwelling longwave flux at surface" ${filenamec}
    #ncatted -a long_name,RH200,m,c,"Relative humidity at 200 mbar pressure surface" -a long_name,RH500,m,c,"Relative humidity at 500 mbar pressure surface" -a long_name,RH700,m,c,"Relative humidity at 700 mbar pressure surface" -a long_name,RH850,m,c,"Relative humidity at 850 mbar pressure surface" -a long_name,OMEGA200,m,c,"Vertical velocity at 200 mbar pressure surface" -a long_name,OMEGA500,m,c,"Vertical velocity at 500 mbar pressure surface" -a long_name,OMEGA700,m,c,"Vertical velocity at 700 mbar pressure surface" -a long_name,OMEGA850,m,c,"Vertical velocity at 850 mbar pressure surface" -a long_name,Z200,m,c,"Geopotential Z at 200 mbar pressure surface" -a long_name,Z500,m,c,"Geopotential Z at 500 mbar pressure surface" -a long_name,Z700,m,c,"Geopotential Z at 700 mbar pressure surface" -a long_name,Z850,m,c,"Geopotential Z at 850 mbar pressure surface" ${filenamec}  
    #ncatted -a long_name,PS,m,c,"Surface pressure" -a long_name,PSL,m,c,"Sea level pressure" -a long_name,TMNUMLIQ,m,c,"NUMLIQ column burden" -a long_name,TMNUMICE,m,c,"NUMICE column burden" -a long_name,TMNUMRAI,m,c,"NUMRAI column burden" ${filenamec}
    #ncatted -a long_name,U,m,c,"Zonal wind" -a long_name,V,m,c,"Meridional wind" ${filenamec}
    #ncatted -a long_name,Q,m,c,"Specific humidity" ${filenamec}
    #ncatted -a long_name,CLDLIQ,m,c,"Grid box averaged cloud liquid amount" -a long_name,CLDICE,m,c,"Grid box averaged cloud ice amount" ${filenamec}
    #ncatted -a long_name,CLOUD,m,c,"Cloud fraction" -a long_name,OMEGA,m,c,"Vertical velocity (pressure)" ${filenamec}
    ncatted -a long_name,EMIS,m,c,"Cloud longwave emissivity" -a long_name,TOT_ICLD_VISTAU,m,c,"Total in-cloud visible optical depth" ${filenamec}
    ncatted -a mixing_ratio,,d,, ${filenamec}
    echo " "
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
