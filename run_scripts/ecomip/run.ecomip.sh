#!/bin/bash
set -e
umask 022

script_root=${PWD}
branch="brhillman" #"decadal-production-run4" #"fix-nanobug-in-horiz-remap" #"decadal-production-run4" #"scorpio-update"
code_root=${HOME}/codes/e3sm/branches/${branch}
res=ne4pg2_ne4pg2 #ne1024pg2_ne1024pg2  #ne1024pg2_ne1024pg2 #ne1024pg2_ne1024pg2 #ne30pg2_EC30to60E2r2
compset=F2010-SCREAMv1
machine="pm-gpu" #chrysalis
compiler="gnugpu" #intel
project="e3sm"
walltime="00:10:00"
datestring="20250929"
casename=${branch}.${res}.${compset}.${datestring}
caseroot=${SCRATCH}/e3sm/cases/${casename}
pecount=32x1

mkdir -p `dirname ${caseroot}`
${code_root}/cime/scripts/create_newcase \
    --case=${caseroot} \
    --res=${res} \
    --compset=${compset} \
    --machine=${machine} \
    --compiler=${compiler} \
    --pecount=${pecount} \
    --project=${project} \
    --walltime=${walltime}

cd ${caseroot}

# Extract input_data_dir for user edits to the namelist
DIN_LOC_ROOT=`./xmlquery DIN_LOC_ROOT --value`

# Change run length
./xmlchange STOP_OPTION=ndays,STOP_N=3
./xmlchange REST_OPTION=ndays,REST_N=3
./xmlchange RESUBMIT=0
./xmlchange RUN_STARTDATE="2024-08-01"

# Run setup before configuring components
./case.setup

# Namelist options for EAMxx
if [ "${res}" == "ne1024pg2_ne1024pg2" ]; then
    ./atmchange initial_conditions::Filename="${DIN_LOC_ROOT}/atm/scream/init/screami_ne1024np4L128_era5-19941001-topoadjx6t_20240214.nc"
elif [ "${res}" == "ne256pg2_ne256pg2" ]; then
    ./atmchange initial_conditions::Filename="${DIN_LOC_ROOT}/atm/scream/init/screami_ne256np4L128_era5-19941001-topoadjx6t_20240123.nc"
fi

# Turn on cosp and set frequency
./atmchange physics::atm_procs_list="mac_aero_mic,rrtmgp,cosp"
./atmchange physics::cosp::cosp_frequency_units="hours"
./atmchange physics::cosp::cosp_frequency=1

output_yaml_files=(`ls ${script_root}/*.yaml`)
for file in ${output_yaml_files[@]}; do
    cp -v ${file} ./
    ./atmchange output_yaml_files+=${file}
done

# Namelist options for ELM
# Specify land initial condition and surface datasets
if [ "${res}" == "ne1024pg2_ne1024pg2" ]; then
cat << EOF >> user_nl_elm
  ! Set input data paths
  finidat = "${DIN_LOC_ROOT}/lnd/clm2/initdata/20231226.I2010CRUELM.ne1024pg2_ICOS10.elm.r.1994-10-01-00000.nc"
EOF
elif [ "${res}" == "ne256pg2_ne256pg2" ]; then
cat << EOF >> user_nl_elm
  ! Set input data paths
  finidat = "${DIN_LOC_ROOT}/lnd/clm2/initdata/20240104.I2010CRUELM.ne256pg2.elm.r.1994-10-01-00000.nc"
EOF
elif [ "${res}" == "ne30pg2_EC30to60E2r2" ]; then
cat << EOF >> user_nl_elm
  finidat = "$DIN_LOC_ROOT/lnd/clm2/initdata/20210802.ICRUELM-1950.ne30pg2_EC30to60E2r2.elm.r.0051-01-01-00000.nc"
EOF
fi

# Set MOSART initial condition
if [ "${res}" == "ne1024pg2_ne1024pg2" ]; then
cat << EOF >> user_nl_mosart
  finidat_rtm = "${DIN_LOC_ROOT}/rof/mosart/initdata/20231226.I2010CRUELM.ne1024pg2_ICOS10.mosart.r.1994-10-01-00000.nc"
EOF
elif [ "${res}" == "ne256pg2_ne256pg2" ]; then
cat << EOF >> user_nl_mosart
  finidat_rtm = "${DIN_LOC_ROOT}/rof/mosart/initdata/20240104.I2010CRUELM.ne256pg2.mosart.r.1994-10-01-00000.nc"
EOF
fi

# Point to new SST forcing
# TODO: need new SST data for 2024
#./xmlchange --file env_run.xml --id SSTICE_DATA_FILENAME --val "${DIN_LOC_ROOT}/atm/cam/sst/sst_ostia_3600x7200_19940930_20151231_c20240125.nc"
#./xmlchange --file env_run.xml --id SSTICE_GRID_FILENAME --val "${DIN_LOC_ROOT}/ocn/docn7/domain.ocn.3600x7200.230522.nc"
#./xmlchange --file env_run.xml --id SSTICE_YEAR_ALIGN --val 1994
#./xmlchange --file env_run.xml --id SSTICE_YEAR_START --val 1994
#./xmlchange --file env_run.xml --id SSTICE_YEAR_END --val 2015

# Link to rundir
ln -s `./xmlquery --value RUNDIR` run

# Copy runscript to case dir
cp ${script_root}/`basename $0` ./

# Setup, build, run
./case.setup
./case.build
./case.submit -a="--job-name=${casename} -t ${walltime} --mail-type=ALL --mail-user=bhillma@sandia.gov"
echo "caseroot = ${caseroot}"
