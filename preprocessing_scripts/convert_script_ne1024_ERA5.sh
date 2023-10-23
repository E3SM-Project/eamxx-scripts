#!/bin/bash
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 01:30:00
#SBATCH -q bigmem
#SBATCH -A e3sm
#SBATCH -o convert_script_ne1024_ERA5-%j.out

module load python
conda activate nco-cmem


simdate=20160801
topo_smooth=x6t
input_initial_condition=/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme/ifs_oper_T1279_2016080100_mod_subset_to_e3sm_ne1024np4_topoadj-x6t_L128.c102721.nc
topo_file=/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne1024np4pg2_x6t-SGH.nc
#/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme/ifs_oper_T1279_2020012000_mod_subset_to_e3sm_ne1024np4_topoadj-x6t_L128.c062121.nc
#/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme/ifs_oper_T1279_2016080100_mod_subset_to_e3sm_ne1024np4_topoadj_L128_v20191116.nc
#/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme/ifs_oper_T1279_2020012000_mod_subset_to_e3sm_ne1024np4_topoadj_L128.nc

#input_initial_condition=${wuyinroot}/ERA5_UVTQ-${simdate}_00_mod_to_e3sm_ne1024np4_topoadj-x6t_L128.nc
#wuyinroot=/global/homes/w/wlin/pe3sm/data/ERA5/ModelLevel/
outputroot=/global/cfs/cdirs/e3sm/bhillma/scream/data/init
#topo_file=/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne1024np4pg2_16xconsistentSGH_20190528_converted.nc
#/global/homes/w/wlin/pe3sm/data/ERA5/ModelLevel/ERA5_UVTQ-20131001_00_mod_to_e3sm_ne1024np4_topoadj-16x_L128.nc
output_initial_condition=${outputroot}/screami_ne1024np4L128_era5-${simdate}-topoadj${topo_smooth}_`date +%Y%m%d`.nc
if [ -e ${input_initial_condition} ]; then
    echo "convert script: ${input_initial_condition} -> ${output_initial_condition}..."
    srun ./convert_cami_to_screami.py \
        ${input_initial_condition} \
        ${topo_file} \
        /global/cfs/cdirs/e3sm/bhillma/scream/data/init/gas_constituents_ne1024np4L128_20220913.nc \
        ${output_initial_condition}
fi
