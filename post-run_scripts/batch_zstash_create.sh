#!/bin/bash

# Load E3SM Unified
source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh

# Lust of experiments to archive with zstash

EXPS=(\
      master.ne1024pg2_r0125_oRRS18to6v3.F2010-SCREAM-HR-DYAMOND2.cori-knl_intel.1536x8x16.DY2_Nov05branch_SHOC_P3_AB_bugfix.20201105-16 \
    )

# Loop over simulations
for EXP in "${EXPS[@]}"
do
    echo === Archiving ${EXP} ====
    cd /global/cscratch1/sd/terai/E3SM_simulations/${EXP}/run
    mkdir -p zstash
    stamp=`date +%Y%m%d`
    time zstash create -v --hpss=zstash_archive/${EXP} . 2>&1 | tee zstash/zstash_create_${stamp}.log
done
