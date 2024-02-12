#!/bin/bash -fe

main() {

do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

readonly MACHINE="frontier-scream-gpu"
readonly CHECKOUT="whannah_2024-autocal-spinup"
readonly BRANCH="master"
readonly COMPILER="crayclang-scream"
readonly DEBUG_COMPILE=FALSE
readonly Q=regular

readonly COMPSET="F2010-SCREAMv1"
readonly RESOLUTION="ne1024pg2_ne1024pg2"

readonly EXPERIMENT='dy1024'
readonly SHARK='nudgec'
readonly jobname="${EXPERIMENT}.${SHARK}.${RESOLUTION}.${COMPSET}.${CHECKOUT}"

readonly SCREAMDOCS_ROOT="/lustre/orion/cli115/proj-shared/noel/wacmy/scream-docs/v1_output/auto-calibration"
readonly CODE_ROOT="/lustre/orion/cli115/proj-shared/noel/wacmy/${CHECKOUT}"
readonly PROJECT="cli115"

githash_eamxx=`git --git-dir ${CODE_ROOT}/.git rev-parse HEAD`
githash_screamdocs=`git --git-dir /lustre/orion/cli115/proj-shared/noel/wacmy/scream-docs/.git rev-parse HEAD`

readonly CASE_NAME=SCREAM.2024-autocal-00.ne1024pg2
readonly CASE_ROOT="/lustre/orion/cli115/proj-shared/noel/e3sm_scratch/${CHECKOUT}/${EXPERIMENT}/${SHARK}/${CASE_NAME}"

readonly CASE_GROUP=""

readonly HIST_OPTION="never"
readonly HIST_N="1"

readonly MODEL_START_TYPE="initial"  # "initial", "continue", "branch", "hybrid"
readonly START_DATE="2020-01-20"     # "" for default, or explicit "0001-01-01"

readonly GET_REFCASE=false
readonly RUN_REFDIR=""
readonly RUN_REFCASE=""
readonly RUN_REFDATE=""   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive
readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
readonly CASE_RUN_DIR=${CASE_ROOT}/run

readonly PELAYOUT="4096x6" # 512 nodes 8 MPI's per node
#readonly PELAYOUT="16384x6" # 2048 nodes 8 MPI's per node
readonly WALLTIME="00:59:00"

readonly STOP_OPTION="ndays"
readonly STOP_N="1"
readonly REST_OPTION="ndays"
readonly REST_N="1"
readonly RESUBMIT="0"
readonly DO_SHORT_TERM_ARCHIVING=false
readonly OLD_EXECUTABLE=""

umask 022

create_newcase

case_setup

case_build

runtime_options

# Copy script into case_script directory for provenance
copy_script

case_submit

echo $'\n----- All done -----\n'

}

user_nl() {

    echo "+++ Configuring SCREAM for 128 vertical levels +++"
    ./xmlchange SCREAM_CMAKE_OPTIONS="SCREAM_NP 4 SCREAM_NUM_VERTICAL_LEV 128 SCREAM_NUM_TRACERS 10"

}

#-----------------------------------------------------
create_newcase() {

    if [ "${do_create_newcase,,}" != "true" ]; then
	echo $'\n----- Skipping create_newcase -----\n'
	return
    fi

    echo $'\n----- Starting create_newcase -----\n'

    args=" --case ${CASE_NAME} \
	--output-root ${CASE_ROOT} \
	--script-root ${CASE_SCRIPTS_DIR} \
	--handle-preexisting-dirs u \
	--compset ${COMPSET} \
	--res ${RESOLUTION} \
	--machine ${MACHINE} \
	--compiler ${COMPILER} \
	--walltime ${WALLTIME} \
	--pecount ${PELAYOUT}"

    if [ ! -z "${PROJECT}" ]; then
      args="${args} --project ${PROJECT}"
    fi
    if [ ! -z "${CASE_GROUP}" ]; then
      args="${args} --case-group ${CASE_GROUP}"
    fi
    if [ ! -z "${QUEUE}" ]; then
      args="${args} --queue ${QUEUE}"
    fi

    ${CODE_ROOT}/cime/scripts/create_newcase ${args}

    if [ $? != 0 ]; then
      echo $'\nNote: if create_newcase failed because sub-directory already exists:'
      echo $'  * delete old case_script sub-directory'
      echo $'  * or set do_newcase=false\n'
      exit 35
    fi

}

#-----------------------------------------------------
case_setup() {

    if [ "${do_case_setup,,}" != "true" ]; then
	echo $'\n----- Skipping case_setup -----\n'
	return
    fi

    echo $'\n----- Starting case_setup -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    ./xmlchange EXEROOT=${CASE_BUILD_DIR}
    ./xmlchange RUNDIR=${CASE_RUN_DIR}
    ./xmlchange DOUT_S=${DO_SHORT_TERM_ARCHIVING}
    ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # Custom user_nl
    user_nl

    ./xmlchange --file env_mach_pes.xml MAX_MPITASKS_PER_NODE="8"
    ./xmlchange --file env_mach_pes.xml MAX_TASKS_PER_NODE="128"
    ./xmlchange --file env_mach_pes.xml NTHRDS="1"
    #./xmlchange --file env_mach_pes.xml NTHRDS_LND="6"
    #./xmlchange --file env_mach_pes.xml NTHRDS_ICE="6"

    ./xmlchange PIO_NETCDF_FORMAT="64bit_data"

    ./case.setup --reset

    # Save provenance info
    echo "branch hash for EAMxx: $githash_eamxx" > GIT_INFO.txt
    echo "master hash for output files: $githash_screamdocs" >> GIT_INFO.txt

    popd
}

#-----------------------------------------------------
case_build() {

    pushd ${CASE_SCRIPTS_DIR}

    # do_case_build = false
    if [ "${do_case_build,,}" != "true" ]; then

	echo $'\n----- case_build -----\n'

	if [ "${OLD_EXECUTABLE}" == "" ]; then
	    # Uses previously built executable, make sure it exists
	    if [ -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
		echo 'Skipping build because $do_case_build = '${do_case_build}
	    else
		echo 'ERROR: $do_case_build = '${do_case_build}' but no executable exists for this case.'
		exit 297
	    fi
	else
	    # If absolute pathname exists and is executable, reuse pre-exiting executable
	    if [ -x ${OLD_EXECUTABLE} ]; then
		echo 'Using $OLD_EXECUTABLE = '${OLD_EXECUTABLE}
		cp -fp ${OLD_EXECUTABLE} ${CASE_BUILD_DIR}/
	    else
		echo 'ERROR: $OLD_EXECUTABLE = '$OLD_EXECUTABLE' does not exist or is not an executable file.'
		exit 297
	    fi
	fi
	echo 'WARNING: Setting BUILD_COMPLETE = TRUE.  This is a little risky, but trusting the user.'
	./xmlchange BUILD_COMPLETE=TRUE

    # do_case_build = true
    else

	echo $'\n----- Starting case_build -----\n'

	if [ "${DEBUG_COMPILE}" == "TRUE" ]; then
	    ./xmlchange DEBUG=${DEBUG_COMPILE}
	fi

	# Run CIME case.build
	./case.build

	# Some user_nl settings won't be updated to *_in files under the run directory
	# Call preview_namelists to make sure *_in and user_nl files are consistent.
	./preview_namelists

    fi

    popd
}

#-----------------------------------------------------
runtime_options() {

    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Set simulation start date
    if [ ! -z "${START_DATE}" ]; then
	./xmlchange RUN_STARTDATE=${START_DATE}
    fi
    # Set temperature cut off in dycore threshold to 180K
    ./atmchange vtheta_thresh=180
    ./atmquery vtheta_thresh

    # Turn on cosp and set default frequency
    ./atmchange physics::atm_procs_list="mac_aero_mic,rrtmgp,cosp,nudging"
    ./atmchange physics::cosp::cosp_frequency_units="hours"
    ./atmchange physics::cosp::cosp_frequency=3
    ./atmchange initial_conditions::Filename="/lustre/orion/cli115/world-shared/e3sm/inputdata/atm/scream/init/screami_ne1024np4L128_ifs-20200120-topoadjx6t_20221011.nc"

    ./atmchange physics::mac_aero_mic::shoc::compute_tendencies=T_mid,qv,qc
    ./atmchange physics::mac_aero_mic::p3::compute_tendencies=T_mid,qv,qc,qr,qi
    ./atmchange physics::rrtmgp::compute_tendencies=T_mid
    ./atmchange homme::compute_tendencies=T_mid,qv
    # use GHG levels more appropriate for 2020
    ./atmchange co2vmr=412.5e-6
    ./atmchange ch4vmr=1877.0e-9
    ./atmchange n2ovmr=332.0e-9
    ./atmchange orbital_year=2020
    # use CO2 the same in land model
    ./xmlchange CCSM_CO2_PPMV=412.5
    #write out DAG
    #./atmchange atmosphere_dag_verbosity_level=5

    #specify land IC file
cat << EOF >> user_nl_elm
 finidat='/lustre/orion/cli115/world-shared/e3sm/inputdata/lnd/clm2/initdata/20221218.F2010-CICE.ne30pg2_ne1024pg2.elm.r.2020-01-20-00000.nc'
 hist_empty_htapes=.true.
EOF

cat << EOF >> user_nl_cpl
 ocn_surface_flux_scheme = 2
EOF

    # things that are different in Walters nudging case:

    ./atmchange lambda_high=0.08

    #nudge_data_root =  '/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data'
    #nudge_map_root  =  '/lustre/orion/cli115/proj-shared/hannah6/HICCUP/files_map'

    ./atmchange nudging::nudging_filename="/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data/HICCUP.nudging_uv_era5.2020-01-20.ne128pg2.L128.nc"
    ./atmchange nudging::nudging_filename+="/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data/HICCUP.nudging_uv_era5.2020-01-21.ne128pg2.L128.nc"
    ./atmchange nudging::nudging_filename+="/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data/HICCUP.nudging_uv_era5.2020-01-22.ne128pg2.L128.nc"
    ./atmchange nudging::nudging_filename+="/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data/HICCUP.nudging_uv_era5.2020-01-23.ne128pg2.L128.nc"
    ./atmchange nudging::nudging_filename+="/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data/HICCUP.nudging_uv_era5.2020-01-24.ne128pg2.L128.nc"
    ./atmchange nudging::nudging_filename+="/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data/HICCUP.nudging_uv_era5.2020-01-25.ne128pg2.L128.nc"
    ./atmchange nudging::nudging_filename+="/lustre/orion/cli115/proj-shared/hannah6/scream_scratch/nudge_data/HICCUP.nudging_uv_era5.2020-01-26.ne128pg2.L128.nc"

     #/lustre/orion/cli115/proj-shared/hannah6/HICCUP/files_map/map_ne128pg2_to_ne1024pg2_trbilin.20240201.nc

    ./atmchange nudging::source_pressure_type=TIME_DEPENDENT_3D_PROFILE
    ./atmchange nudging::nudging_fields=U,V
    ./atmchange nudging::nudging_timescale=21600 #  6-hr
    ./atmchange nudging::nudging_refine_remap_mapfile='/lustre/orion/cli115/proj-shared/hannah6/HICCUP/files_map/map_ne128pg2_to_ne1024pg2_trbilin.20240201.nc'

    # Segment length
    ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}

    # Restart frequency
    ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}

    # Coupler history
    ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}

    # Coupler budgets (always on)
    ./xmlchange BUDGETS=TRUE

    # Set resubmissions
    if (( RESUBMIT > 0 )); then
	./xmlchange RESUBMIT=${RESUBMIT}
    fi

    # Start from default of user-specified initial conditions
    if [ "${MODEL_START_TYPE,,}" == "initial" ]; then
	./xmlchange RUN_TYPE="startup"
	./xmlchange CONTINUE_RUN="FALSE"

    # Continue existing run
    elif [ "${MODEL_START_TYPE,,}" == "continue" ]; then
	./xmlchange CONTINUE_RUN="TRUE"

    elif [ "${MODEL_START_TYPE,,}" == "branch" ] || [ "${MODEL_START_TYPE,,}" == "hybrid" ]; then
	./xmlchange RUN_TYPE=${MODEL_START_TYPE,,}
	./xmlchange GET_REFCASE=${GET_REFCASE}
	./xmlchange RUN_REFDIR=${RUN_REFDIR}
	./xmlchange RUN_REFCASE=${RUN_REFCASE}
	./xmlchange RUN_REFDATE=${RUN_REFDATE}
	echo 'Warning: $MODEL_START_TYPE = '${MODEL_START_TYPE}
	echo '$RUN_REFDIR = '${RUN_REFDIR}
	echo '$RUN_REFCASE = '${RUN_REFCASE}
	echo '$RUN_REFDATE = '${START_DATE}

    else
	echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
	exit 380
    fi



    ./xmlchange --file env_run.xml --id SSTICE_DATA_FILENAME --val "/lustre/orion/cli115/world-shared/e3sm/inputdata/atm/cam/sst/sst_ifs_360x720_20200110_20200305_c201013.nc"
    ./xmlchange --file env_run.xml --id SSTICE_GRID_FILENAME --val "/lustre/orion/cli115/world-shared/e3sm/inputdata/ocn/docn7/domain.ocn.360x720.201027.nc"
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_ALIGN --val 2020
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_START --val 2020
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_END --val 2020


    popd
}

#-----------------------------------------------------
case_submit() {

    if [ "${do_case_submit,,}" != "true" ]; then
	echo $'\n----- Skipping case_submit -----\n'
	return
    fi

    echo $'\n----- Starting case_submit -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Run CIME case.submit
    ./case.submit -a="--job-name=${jobname} -t ${WALLTIME} --mail-type=END --mail-user=ndkeen@lbl.gov" >& submitout.txt
    #./case.submit -a="--qos=${Q}"

    popd
}

#-----------------------------------------------------
copy_script() {

    echo $'\n----- Saving run script for provenance -----\n'

    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=`basename $0`
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp ${this_script_name} ${script_provenance_dir}/${script_provenance_name}

}

#-----------------------------------------------------
# Silent versions of popd and pushd
pushd() {
    command pushd "$@" > /dev/null
}
popd() {
    command popd "$@" > /dev/null
}

# Now, actually run the script
#-----------------------------------------------------
main
