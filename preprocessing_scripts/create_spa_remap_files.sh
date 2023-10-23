#!/bin/bash
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 00:30:00
#SBATCH -q bigmem
#SBATCH -A e3sm
#SBATCH -o create_spa_remap_files-%j.out

datestr=20221011
algorithm="intbilin_se2fv"

mapping_root=/global/cfs/cdirs/e3sm/mapping
output_root=/global/cfs/cdirs/e3sm/inputdata/atm/scream/maps

source_grid_name=ne30np4
source_grid_file=${mapping_root}/grids/ne30.g

target_grid_name=ne1024pg2
target_grid_file=${mapping_root}/grids/ne1024pg2_scrip_20221011.nc

map_file=${output_root}/map_${source_grid_name}_to_${target_grid_name}_${algorithm}_${datestr}.nc
ncremap -5 -a ${algorithm} -s ${source_grid_file} -g ${target_grid_file} -m ${map_file}
