#!/bin/bash
set -e
umask 022

script_root=${PWD}
branch="decadal-production-run6" #"decadal-production-run4" #"fix-nanobug-in-horiz-remap" #"decadal-production-run4" #"scorpio-update"
code_root=${HOME}/codes/scream/branches/${branch}
screamdocs_root=${HOME}/codes/scream-decadal/scream-docs
res=ne1024pg2_ne1024pg2  #ne1024pg2_ne1024pg2 #ne1024pg2_ne1024pg2 #ne30pg2_EC30to60E2r2
compset=F20TR-SCREAMv1
machine="frontier-scream-gpu" #chrysalis
compiler="craygnuamdgpu" #"crayclang-scream" #"craygnuamdgpu" #"crayclang-scream" #intel
project="cli115"
walltime="01:00:00"
datestring="20240708"
RUN_REFDATE="1999-04-05"; RUN_REFTOD="00000"
#RUN_REFDATE="1998-09-01"; RUN_REFTOD="00000"
casename=${branch}-${datestring}.${res}.${compset}.pnetcdf # NOTE: cannot change this for pseudo-branch
casebase=${compiler}-branch # Need a unique id to tell us where to write run and bld so we do not overwrite original
caseroot=${HOME}/codes/scream/cases/scream-decadal/${casebase}/${casename}
#readonly pecount="1536x6" # 192 nodes
#readonly pecount="3072x6" # 384 nodes
#readonly pecount="4096x6" # 512 nodes
#readonly pecount="8192x6" # 1024 nodes
#readonly pecount="15056x6" # 1882 nodes
if [ "${res}" == "ne1024pg2_ne1024pg2" ]; then
    #pecount="2560x6" # 320 nodes 
    pecount="16384x6" # 2048 nodes
elif [ "${res}" == "ne256pg2_ne256pg2" ]; then
    pecount="768x6" # 96 nodes
elif [ "${res}" == "ne30pg2_EC30to60E2r2" ]; then
    pecount="16x6"
    #pecount=5540x1 #128x1
fi
mkdir -p `dirname ${caseroot}`
if [ 1 -eq 1 ]; then
${code_root}/cime/scripts/create_newcase \
    --case=${caseroot} \
    --res=${res} \
    --compset=${compset} \
    --machine=${machine} \
    --compiler=${compiler} \
    --pecount=${pecount} \
    --project=${project} \
    --walltime=${walltime}
fi
cd ${caseroot}

# Extract input_data_dir for user edits to the namelist
DIN_LOC_ROOT=`./xmlquery DIN_LOC_ROOT --value`

# Change run and bld directories so we can use the same casename for a pseudo-branch run
./xmlchange RUNDIR=/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/scream-decadal/${casebase}/${casename}/run
./xmlchange EXEROOT=/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/scream-decadal/${casebase}/${casename}/bld

# Link to rundir
#ln -s `./xmlquery --value RUNDIR` run

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
    #./case.setup --reset
fi

# Change run length
./xmlchange STOP_OPTION=ndays,STOP_N=1
./xmlchange REST_OPTION=ndays,REST_N=1
./xmlchange RESUBMIT=0
./xmlchange RUN_STARTDATE="1994-10-01"

# Configure for (pseudo) branch run
if [ 1 -eq 1 ]; then
    RUN_REFDIR="/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-run6-20240708.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.pnetcdf/run"
    ./xmlchange RUN_TYPE="startup"
    #./xmlchange RUN_REFCASE="decadal-production-run6-20240708.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.pnetcdf"
    #./xmlchange RUN_REFDIR=${RUN_REFDIR}
    #./xmlchange RUN_REFDATE=${RUN_REFDATE}
    #./xmlchange RUN_REFTOD=${RUN_REFTOD}
    #./xmlchange GET_REFCASE="FALSE"
    ./xmlchange CONTINUE_RUN=TRUE 
    # Copy files ourselves
    RUNDIR=`./xmlquery --value RUNDIR`
    mkdir -p ${RUNDIR}
    cp ${RUN_REFDIR}/rpointer.??? ${RUNDIR}
    for file in ${RUNDIR}/rpointer.???; do
        sed -i "s/[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]-00000/${RUN_REFDATE}-${RUN_REFTOD}/" ${file}
        for f in `cat ${file}`; do
            echo ${f}
            if [ -f ${RUN_REFDIR}/${f} ] && [ ! -f ${RUNDIR}/`basename ${f}` ]; then
                ln -s ${RUN_REFDIR}/`basename ${f}` ${RUNDIR}/`basename ${f}`
            fi
        done
    done
    # Also need to copy the rh files, which might not be picked up in rpointer files
    for file in ${RUN_REFDIR}/*.rh[0-9].${RUN_REFDATE}-${RUN_REFTOD}.nc; do
        if [ ! -f ${RUNDIR}/`basename ${file}` ]; then
            ln -s ${file} ${RUNDIR}/
        fi
    done
    # And copy .bin file from correct dummy case
    RUN_REFDIR="/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/scream-decadal/craygnuamdgpu-gen-docn/decadal-production-run6-20240708.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.pnetcdf/run"
    for file in ${RUN_REFDIR}/*.${RUN_REFDATE}-${RUN_REFTOD}.bin; do
        if [ -e ${RUNDIR}/`basename ${file}` ]; then
            rm ${RUNDIR}/`basename ${file}`
        fi
        ln -s ${file} ${RUNDIR}/
    done
fi

# Turn on budget reporting
./xmlchange BUDGETS=TRUE

# For big data
./xmlchange PIO_NETCDF_FORMAT="64bit_data"
./xmlchange PIO_TYPENAME=pnetcdf #adios #,PIO_TYPENAME_ATM=adios
./xmlchange PIO_REARRANGER=1  # use PIO_REARRANGER=3, for ADIOS; PIO_REARRANGER=1 for pnetcdf
if [ "${pecount}" == "2560x6" ]; then
    ./xmlchange PIO_STRIDE=4  #BRHDEBUG only needed for small node counts
fi

# Change to non-parallel netcdf
#sed -i 's|<command name="load">cray-hdf5-parallel.*|<command name="load">cray-hdf5/1.12.1.5</command>|' env_mach_specific.xml
#sed -i 's|<command name="load">cray-netcdf-hdf5parallel.*|<command name="load">cray-netcdf/4.8.1.5</command>|' env_mach_specific.xml

# Run setup before configuring components
./case.setup

# Automatic restarts on regex
#./xmlchange NODE_FAIL_REGEX=''
#./xmlchange ALLOCATE_SPARE_NODES=TRUE
#./xmlchange MPIRUN_RETRY_REGEX=''
#./xmlchange MPIRUN_RETRY_COUNT=1

# Extra diagnostics for non-determinism debugging
#./atmchange --all \
#    internal_diagnostics_level=1 \
#    atmosphere_processes::internal_diagnostics_level=0 \
#    ctl_nl::internal_diagnostics_level=0

# Namelist options for EAMxx
if [ "${res}" == "ne1024pg2_ne1024pg2" ]; then
    ./atmchange initial_conditions::Filename="${DIN_LOC_ROOT}/atm/scream/init/screami_ne1024np4L128_era5-19941001-topoadjx6t_20240214.nc"
elif [ "${res}" == "ne256pg2_ne256pg2" ]; then
    ./atmchange initial_conditions::Filename="${DIN_LOC_ROOT}/atm/scream/init/screami_ne256np4L128_era5-19941001-topoadjx6t_20240123.nc"
fi

# Run with bugfixed SPA file
./atmchange spa_data_file="${DIN_LOC_ROOT}/atm/scream/init/spa_v3.LR.F2010.2011-2025.c_20240405.nc"

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
    sed -i "s|horiz_remap_file:.*_to_ne30.*|horiz_remap_file: ${map_to_ne30}|" ./`basename ${file}`
    sed -i "s|horiz_remap_file:.*_to_DecadalSites.*|horiz_remap_file: ${map_to_DecadalSites}|" ./`basename ${file}`
done

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
  hist_mfilt = 1,1
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

# Copy runscript to case dir
cp ${script_root}/`basename $0` ./

# Setup, build, run
./case.setup
./case.build
#./case.submit -a="--job-name=decadal-amip -t 12:00:00  --mail-type=ALL --mail-user=bhillma@sandia.gov -x frontier08656 --account=cli115" # --gpu-srange=800-1600"
./case.submit -a="--job-name=decadal-amip -t ${walltime} --qos=debug --mail-type=ALL --mail-user=bhillma@sandia.gov -x frontier08656 --account=cli115" # --gpu-srange=800-1600"
echo "${caseroot}"
