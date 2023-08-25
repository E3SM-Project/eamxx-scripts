#!/bin/bash                                 
#$ -cwd                                     

# regrids ne1024np4 output to 70 km fv grid 
# Tested script having requested 1 interactive nodes on cori large-memory nodes (amd)
# Instructions courtesy of Charlie Zender
# Mainly for regridding 2D variables        

# ============================================================                                                

# Location of the model output              
# USER needs to change this                 
fileloc='/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127/'

# Map file
map_file='/global/cfs/cdirs/e3sm/mapping/maps/map_ne1024pg2_to_fv256x512_nco.20201201.nc'

# Alternative map file for high resolution regridding
#map_file=/global/cfs/cdirs/e3sm/terai/mapping/'map_ne1024pg2_to_fv1800x3600_nco.20210407.nc'


# Location of the regridded files
# USER needs to change this
outputloc='/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127/regridded/'


# =============================================================                                               

# Calls a bash script that retrieves, regrids, and concatenates each variable data                            
# Variables can be added or removed based on the interest of the USER, but currently regrids h0, h1, h2, h3, h4 files in order                                                  

# USER: Compare the fincl list with the list below to make sure they correspond
echo "  "
echo " ---------------- Start regridding -------------------- "
echo "  "


cd ${fileloc}

ls ${casename}*eam.h?.2020-01-27-00000.nc | ncremap --dbg=1 --vrb=3 --devnull=Np --nco='--dbg=5' --thr_nbr=3 --par_typ=bck --job_nbr=10 -m ${map_file} -O ${outputloc} > ./ncremap.2020-01-27 2>&1 
#echo "done with 02-24"


echo "  "
echo "  "
echo "  "
echo "   DONE!  "
echo "  "
echo "  "
