#!/bin/bash

#SBATCH --account=e3sm
#SBATCH --job-name moab
#SBATCH --constraint=cpu
#SBATCH --time=00:30:00
#SBATCH -q debug     
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --nodes=8     
#SBATCH --ntasks-per-node=128     
#SBATCH --cpus-per-task=1       


source /global/u1/m/meng/utils/e3sm_mach/pm-cpu_intel.sh 
moabexe=/global/cfs/cdirs/e3sm/software/moab/intel-develop/bin
mybin=/global/homes/m/meng/.conda/envs/my_root/bin

n=1024
res='1024x4'
mapfile="map_721x1440_to_conus26_ne${res}np4.nc"

dst="/global/cfs/cdirs/e3sm/2026-INCITE-CONUS-RRM/files_grid/2026-incite-conus-${res}.g"
src="/pscratch/sd/m/meng/hiccup/incite_conus/scrip_ERA5_721x1440.nc"

dst_h5m="/pscratch/sd/m/meng/hiccup/incite_conus/tmp/2026-incite-conus-${res}"
src_h5m="/pscratch/sd/m/meng/hiccup/incite_conus/tmp/scrip_ERA5_721x1440"

dst_h5m_p=${dst_h5m}_${n}p.h5m
src_h5m_p=${src_h5m}_${n}p.h5m

${moabexe}/mbconvert -B -o PARALLEL=WRITE_PART -O PARALLEL=BCAST_DELETE \
    -O PARTITION=TRIVIAL -O PARALLEL_RESOLVE_SHARED_ENTS \
    $src ${src_h5m}.h5m

${moabexe}/mbconvert -B -o PARALLEL=WRITE_PART -O PARALLEL=BCAST_DELETE \
    -O PARTITION=TRIVIAL -O PARALLEL_RESOLVE_SHARED_ENTS \
    -i GLOBAL_DOFS -r 4 \
    $dst ${dst_h5m}.h5m

${moabexe}/mbpart $n --zoltan RCB --globalIds ${src_h5m}.h5m ${src_h5m_p}
${moabexe}/mbpart $n --zoltan RCB --recompute_rcb_box --scale_sphere --project_on_sphere 2 ${dst_h5m}.h5m ${dst_h5m_p}

srun ${moabexe}/mbtempest --type 5 --weights \
    --load ${dst_h5m_p} --load ${src_h5m_p} \
    --method cgll --order 4 --global_id GLOBAL_DOFS --method fv --order 1 --global_id GLOBAL_ID \
    --monotonicity 1 --boxeps 1e-13 --gnomonic --verbose \
    --file /pscratch/sd/m/meng/hiccup/incite_conus/tmp/map_tmp${res}.nc --sparseconstraints

GenerateTransposeMap --in /pscratch/sd/m/meng/hiccup/incite_conus/tmp/map_tmp${res}.nc \
    --out /pscratch/sd/m/meng/hiccup/incite_conus/$mapfile 

exit 0 