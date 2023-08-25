#!/bin/bash
#$ -cwd

# regrids ne1024np4 output to 70 km fv grid
# Tested script having requested 4 interactive nodes on cori-haswell
# This script has a different formulation than the ne120 regridding script
# because each output file is too big.

# We break up the process by regridding specific variables at a time from h1/h2 file.
# {1} = input directory with se model output data
# {2} = variable name
# {3} = htype (e.g., h0, h1, h2, ...)
# {4} = output directory for regridded data

# Note: this will likely fail with memory issues in interactive nodes
# if 3D variables with multiple timesteps are processed. 


echo " ${1} "

cd ${1}

echo "Regridding ${2} variable from ${3} files"

for f in *.eam.${3}.*.nc
do
    echo "processing $f"
    ncremap -v ${2} -m /global/cfs/cdirs/e3sm/terai/mapping/map_ne1024pg2_to_fv256x512_nco.20201201.nc -i ${1}${f} -o ${4}${2}_rgr_256x512_${f}
done

ncrcat ${4}${2}_rgr_256x512_*nc ${4}CAT_${2}_rgr_256x512_${3}.nc

echo "  "
echo "  "
echo "  "
echo "   DONE!  "
echo "  "
echo "  "
echo "  ******** "  
