#!/bin/bash
set -e
umask 022

script_root=${PWD}
branch="add-F20TR-SCREAMv1"
code_root=${HOME}/codes/scream/branches/${branch}
screamdocs_root=${HOME}/codes/scream-docs/branches/add-decadal-outputs
res=ne256pg2_ne256pg2 #ne1024pg2_ne1024pg2 #ne1024pg2_ne1024pg2 #ne30pg2_EC30to60E2r2
compset=F20TR-SCREAMv1
machine="frontier-scream-gpu" #chrysalis
compiler="crayclang-scream" #intel
project="cli115"
walltime="01:00:00"
datestring="20240216"
casename=decadal-amip.${res}.${compset}.${datestring}-testRestartFixed
caseroot=${HOME}/codes/scream/cases/${casename}
#readonly pecount="1536x6" # 192 nodes
#readonly pecount="3072x6" # 384 nodes
#readonly pecount="4096x6" # 512 nodes
#readonly pecount="8192x6" # 1024 nodes
#readonly pecount="15056x6" # 1882 nodes
if [ "${res}" == "ne1024pg2_ne1024pg2" ]; then
    pecount="2560x6" # 320 nodes 
    #pecount="16384x6" # 2048 nodes
elif [ "${res}" == "ne256pg2_ne256pg2" ]; then
    pecount="768x6" # 96 nodes
elif [ "${res}" == "ne30pg2_EC30to60E2r2" ]; then
    pecount="16x6"
    #pecount=5540x1 #128x1
fi
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

# Change threads
if [ "${machine}" == "frontier-scream-gpu" ]; then
    ./xmlchange --file env_mach_pes.xml NTHRDS="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ATM="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_LND="6"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ICE="6"
    ./xmlchange --file env_mach_pes.xml NTHRDS_OCN="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ROF="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_CPL="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_GLC="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_WAV="1"
fi

# Change run length
./xmlchange STOP_OPTION=nhours,STOP_N=1
./xmlchange REST_OPTION=nhours,REST_N=1
./xmlchange RESUBMIT=1
./xmlchange RUN_STARTDATE="1994-10-01"

# For big data
./xmlchange PIO_NETCDF_FORMAT="64bit_data"
if [ "${pecount}" == "2560x6" ]; then
    ./xmlchange PIO_STRIDE=4  #BRHDEBUG only needed for small node counts
fi

# Run setup before configuring components
./case.setup

# Automatic restarts on regex
#./xmlchange NODE_FAIL_REGEX=''
#./xmlchange ALLOCATE_SPARE_NODES=TRUE
#./xmlchange MPIRUN_RETRY_REGEX=''
#./xmlchange MPIRUN_RETRY_COUNT=1

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

# Need to explicitly turn on computing tendencies
./atmchange physics::mac_aero_mic::shoc::compute_tendencies=T_mid,qv
./atmchange physics::mac_aero_mic::p3::compute_tendencies=T_mid,qv
./atmchange physics::rrtmgp::compute_tendencies=T_mid
./atmchange homme::compute_tendencies=T_mid,qv

# Set temperature cut off in dycore threshold to 180K
./atmchange vtheta_thresh=180

# Change lambda_high
./atmchange lambda_high=0.08

# use GHG levels more appropriate for sim
# Average from 19940101 - 20150101
./atmchange co2vmr=377.2e-6
./atmchange ch4vmr=1786.6e-9
./atmchange n2ovmr=318.6e-9
./atmchange orbital_year=-9999
# use CO2 the same in land model
./xmlchange CCSM_CO2_PPMV=377.2
# write out DAG
./atmchange atmosphere_dag_verbosity_level=5

# Copy output stream yaml files
if [ "${res}" == "ne1024pg2_ne1024pg2" ]; then
    map_to_ne30="${DIN_LOC_ROOT}/atm/scream/maps/map_ne1024pg2_to_ne30pg2_mono.20230901.nc"
    map_to_DecadalSites="${DIN_LOC_ROOT}/atm/scream/maps/map_ne1024pg2_to_DecadalSites_c20240130.nc"
elif [ "${res}" == "ne256pg2_ne256pg2" ]; then
    map_to_ne30="${DIN_LOC_ROOT}/atm/scream/maps/map_ne256pg2_to_ne30pg2_traave.20240206.nc"
    map_to_DecadalSites="${DIN_LOC_ROOT}/atm/scream/maps/map_ne256pg2_to_DecadalSites_c20240130.nc"
else
    echo "Unsupported res for horiz maps"
    exit 1
fi
output_yaml_files=(`ls ${screamdocs_root}/v1_output/decadal/*.yaml`)
for file in ${output_yaml_files[@]}; do
    # TODO: add remap file replacement for different grids
    cp -v ${file} ./
    if [ "${file}" == "${output_yaml_files[0]}" ]; then
        # First file, reset output list
        ./atmchange output_yaml_files="./`basename ${file}`"
    else
        # Append to output list
        ./atmchange output_yaml_files+="./`basename ${file}`"
    fi
    # Replace remap files
    sed -i "s|horiz_remap_file:*._to_ne30.*|horiz_remap_file: ${map_to_ne30}|" ./`basename ${file}`
    sed -i "s|horiz_remap_file:*._to_DecadalSites.*|horiz_remap_file: ${map_to_DecadalSites}|" ./`basename ${file}`
done
#BRHDEBUG
#./atmchange output_yaml_files="./output_debug.yaml"
#cat << EOF >> ./output_debug.yaml
#%YAML 1.1
#---
#filename_prefix: output.scream.decadal.monthlyAVG_ne30pg2
#Averaging Type: Average
#Max Snapshots Per File: 1
#horiz_remap_file: ${DIN_LOC_ROOT}/atm/scream/maps/map_ne1024pg2_to_ne30pg2_mono.20230901.nc
#Fields:
#  Physics PG2:
#    Field Names:
#    - PotentialTemperature_at_700hPa
#output_control:
#  Frequency: 1
#  frequency_units: nmonths
#  MPI Ranks in Filename: false
#EOF
#BRHDEBUG

# Namelist options for ELM
# Specify land initial condition and surface datasets
if [ "${res}" == "ne1024pg2_ne1024pg2" ]; then
cat << EOF >> user_nl_elm
  ! Set input data paths
  finidat = "${DIN_LOC_ROOT}/lnd/clm2/initdata/20231226.I2010CRUELM.ne1024pg2_ICOS10.elm.r.1994-10-01-00000.nc"
  flanduse_timeseries = "${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_ne1024pg2_historical_simyr1990-2014_c240109.nc"
  fsurdat = "${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/surfdata_ne1024pg2_simyr2010_c211021.nc"
EOF
elif [ "${res}" == "ne256pg2_ne256pg2" ]; then
cat << EOF >> user_nl_elm
  ! Set input data paths
  finidat = "${DIN_LOC_ROOT}/lnd/clm2/initdata/20240104.I2010CRUELM.ne256pg2.elm.r.1994-10-01-00000.nc"
  flanduse_timeseries = "${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_ne256pg2_hist_simyr1850-2015_c240131.nc"
  fsurdat = "${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/surfdata_ne256pg2_simyr1850_c240131.nc"
EOF
elif [ "${res}" == "ne30pg2_EC30to60E2r2" ]; then
cat << EOF >> user_nl_elm
  finidat = "$DIN_LOC_ROOT/lnd/clm2/initdata/20210802.ICRUELM-1950.ne30pg2_EC30to60E2r2.elm.r.0051-01-01-00000.nc"
  flanduse_timeseries = "${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_ne30np4.pg2_hist_simyr1850-2015_c210113.nc"
  fsurdat = "${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/surfdata_ne30pg2_simyr1850_c230417_with_TOP.nc"
EOF
fi

# Additional settings for land
# TODO: revisit whether or not we need daily means?
cat << EOF >> user_nl_elm
  ! Override consistency checks since we know our surface data is inconsistent
  check_finidat_fsurdat_consistency = .false.
  check_finidat_year_consistency = .false.
  check_finidat_pct_consistency = .false.
  check_dynpft_consistency = .false.

  ! Make sure we do transient PFTs
  do_transient_pfts = .true.
 
  hist_dov2xy = .true.,.true.
  hist_mfilt = 1,120
  hist_nhtfrq = 0,-24
  hist_avgflag_pertape = 'A','A'
  hist_fincl1 = 'FIRE', 'FPSN', 'QDRAI', 'QRUNOFF', 'ZWT', 'FSAT', 'H2OSOI', 'EFLX_LH_TOT',
                'QVEGT', 'QVEGE', 'FSH', 'ALBD', 'ALBI', 'TBOT', 'QBOT', 'RAIN', 'SNOW',
                'FSDS', 'FSDSND', 'FSDSNI', 'FSDSVD', 'FSDSVI', 'FLDS'
  hist_fincl2 = 'H2OSNO','SOILWATER_10CM','TG'
EOF

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

# Coupler settings; new surface flux scheme
cat << EOF >> user_nl_cpl
  ocn_surface_flux_scheme = 2
EOF

# Point to new SST forcing
./xmlchange --file env_run.xml --id SSTICE_DATA_FILENAME --val "${DIN_LOC_ROOT}/atm/cam/sst/sst_ostia_3600x7200_19940930_20151231_c20240125.nc"
./xmlchange --file env_run.xml --id SSTICE_GRID_FILENAME --val "${DIN_LOC_ROOT}/ocn/docn7/domain.ocn.3600x7200.230522.nc"
./xmlchange --file env_run.xml --id SSTICE_YEAR_ALIGN --val 1994
./xmlchange --file env_run.xml --id SSTICE_YEAR_START --val 1994
./xmlchange --file env_run.xml --id SSTICE_YEAR_END --val 2015

# Link to rundir
ln -s `./xmlquery --value RUNDIR` run

# Copy runscript to case dir
cp ${script_root}/`basename $0` ./

# Setup, build, run
./case.setup
./case.build
./case.submit -a="--job-name=${casename} -t ${walltime} --mail-type=ALL --mail-user=bhillma@sandia.gov"
echo "${caseroot}"
