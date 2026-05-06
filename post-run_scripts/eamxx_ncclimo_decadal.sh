#!/bin/bash                                 

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

drc_in=/pscratch/sd/t/terai/e3sm_scratch/pm-gpu/ne256pg2_ne256pg2.F20TR-SCREAMv1.copilot-fix.cosp-test/monthly/
drc_out=/pscratch/sd/t/terai/e3sm_scratch/pm-gpu/ne256pg2_ne256pg2.F20TR-SCREAMv1.copilot-fix.cosp-test/climo/
caseid=1ma.AVERAGE.nmonths_x1.  #Change this to match the output file prefix

# create climatology files
cd ${drc_in};ls ${caseid}*.nc | ncclimo -P eamxx --fml_nm=eamxx_ne256 --yr_srt=1995 --yr_end=1995 --drc_out=$drc_out --jobs=4

map=/global/cfs/projectdirs/e3sm/zender/maps/map_ne30pg2_to_cmip6_180x360_traave.20231201.nc
#map=/global/cfs/projectdirs/e3sm/zender/maps/map_ne256pg2_to_cmip6_180x360_traave.20250301.nc
cd $drc_out;ls *.nc | ncremap -P eamxx --prm_opt=time,lwband,swband,ilev,lev,plev,cosp_tau,cosp_cth,cosp_prs,dim2,ncol --map=${map} --drc_out=${drc_out}/rgr


