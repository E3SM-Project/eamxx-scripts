#!/bin/bash -fe

# EAMxx template run script

main() {

do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

readonly MACHINE="frontier-scream-gpu"
readonly CHECKOUT="20230916"
readonly BRANCH="simulations/cess-production"
readonly CHERRY=( )
readonly COMPILER="crayclang-scream"
readonly DEBUG_COMPILE=FALSE
readonly Q=regular

# Simulation
readonly COMPSET="F2010-SCREAMv1"
readonly RESOLUTION="ne1024pg2_ne1024pg2"

readonly SCREAMDOCS_ROOT="/ccs/home/terai/scream-docs"
readonly CODE_ROOT="/ccs/home/terai/SCREAM/code_for_frontier"
readonly PROJECT="cli115"

githash_eamxx=`git --git-dir ${CODE_ROOT}/.git rev-parse HEAD`
githash_screamdocs=`git --git-dir ${SCREAMDOCS_ROOT}/.git rev-parse HEAD`

readonly CASE_NAME=cess-cntl.${RESOLUTION}.${COMPSET}.${CHECKOUT}

readonly CASE_ROOT="/lustre/orion/cli115/proj-shared/terai/e3sm_scratch/${CASE_NAME}"

readonly CASE_GROUP=""

# History file frequency (if using default above)
readonly HIST_OPTION="nmonths"
readonly HIST_N="1"

# Run options
readonly MODEL_START_TYPE="initial"  # "initial", "continue", "branch", "hybrid"
readonly START_DATE="2019-08-01"     # "" for default, or explicit "0001-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=false
readonly RUN_REFDIR=""
readonly RUN_REFCASE=""
readonly RUN_REFDATE=""   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'


# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
readonly CASE_RUN_DIR=${CASE_ROOT}/run
#readonly PELAYOUT="1536x6" # 192 nodes
#readonly PELAYOUT="3072x6" # 384 nodes
#readonly PELAYOUT="4096x6" # 512 nodes
#readonly PELAYOUT="8192x6" # 1024 nodes
#readonly PELAYOUT="15056x6" # 1882 nodes
readonly PELAYOUT="16384x6" # 2048 nodes
readonly WALLTIME="00:29:00"
readonly STOP_OPTION="ndays"
readonly STOP_N="1"
readonly REST_OPTION="ndays"
readonly REST_N="1"
readonly RESUBMIT="0"
readonly DO_SHORT_TERM_ARCHIVING=false

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

# Create case
create_newcase

# Setup
case_setup

# Build
case_build

# Configure runtime options
runtime_options

# Copy script into case_script directory for provenance
copy_script

# Submit
#case_submit -a="--qos=${Q}"
case_submit

# All done
echo $'\n----- All done -----\n'

}

# =======================
# Custom user_nl settings
# =======================

user_nl() {

    echo "+++ Configuring SCREAM for 128 vertical levels +++"
    ./xmlchange SCREAM_CMAKE_OPTIONS="SCREAM_NP 4 SCREAM_NUM_VERTICAL_LEV 128 SCREAM_NUM_TRACERS 10"

}

######################################################
### Most users won't need to change anything below ###
######################################################

#-----------------------------------------------------
fetch_code() {

    if [ "${do_fetch_code,,}" != "true" ]; then
	echo $'\n----- Skipping fetch_code -----\n'
	return
    fi

    echo $'\n----- Starting fetch_code -----\n'
    local path=${CODE_ROOT}
    local repo=scream

    echo "Cloning $repo repository branch $BRANCH under $path"
    if [ -d "${path}" ]; then
	echo "ERROR: Directory already exists. Not overwriting"
	exit 20
    fi
    mkdir -p ${path}
    pushd ${path}

    # This will put repository, with all code
    git clone git@github.com:E3SM-Project/${repo}.git .

    # Q: DO WE NEED THIS FOR EAMXX?
    # Setup git hooks
    rm -rf .git/hooks
    git clone git@github.com:E3SM-Project/E3SM-Hooks.git .git/hooks
    git config commit.template .git/hooks/commit.template

    # Check out desired branch
    git checkout ${BRANCH}

    # Custom addition
    if [ "${CHERRY}" != "" ]; then
	echo ----- WARNING: adding git cherry-pick -----
	for commit in "${CHERRY[@]}"
	do
	    echo ${commit}
	    git cherry-pick ${commit}
	done
	echo -------------------------------------------
    fi

    # Bring in all submodule components
    git submodule update --init --recursive

    popd
}

#-----------------------------------------------------
create_newcase() {

    if [ "${do_create_newcase,,}" != "true" ]; then
	echo $'\n----- Skipping create_newcase -----\n'
	return
    fi

    echo $'\n----- Starting create_newcase -----\n'

    # Base arguments
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

    # Oprional arguments
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

    # Setup some CIME directories
    ./xmlchange EXEROOT=${CASE_BUILD_DIR}
    ./xmlchange RUNDIR=${CASE_RUN_DIR}

    # Short term archiving
    ./xmlchange DOUT_S=${DO_SHORT_TERM_ARCHIVING}
    ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # Custom user_nl
    user_nl

    ./xmlchange --file env_mach_pes.xml NTHRDS="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ATM="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_LND="6"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ICE="6"
    ./xmlchange --file env_mach_pes.xml NTHRDS_OCN="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ROF="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_CPL="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_GLC="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_WAV="1"

    ./xmlchange PIO_NETCDF_FORMAT="64bit_data"

    # Finally, run CIME case.setup
    ./case.setup --reset

    # Save provenance invfo
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
	    # Ues previously built executable, make sure it exists
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

	# Turn on debug compilation option if requested
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
    ./atmchange physics::atm_procs_list="mac_aero_mic,rrtmgp,cosp"
    ./atmchange physics::cosp::cosp_frequency_units="hours"
    ./atmchange physics::cosp::cosp_frequency=1

    # Turn on turbulent mountain stress
    ./atmchange physics::mac_aero_mic::atm_procs_list="tms,shoc,cldFraction,spa,p3"

    ./atmchange initial_conditions::Filename="/lustre/orion/cli115/world-shared/e3sm/inputdata/atm/scream/init/screami_ne1024np4L128_era5-20190801-topoadjx6t_20230620.nc"

    ./atmchange physics::mac_aero_mic::shoc::compute_tendencies=T_mid,qv
    ./atmchange physics::mac_aero_mic::p3::compute_tendencies=T_mid,qv
    ./atmchange physics::rrtmgp::compute_tendencies=T_mid
    ./atmchange homme::compute_tendencies=T_mid,qv
    # use GHG levels more appropriate for 2019
    ./atmchange co2vmr=410.5e-6
    ./atmchange ch4vmr=1877.0e-9
    ./atmchange n2ovmr=332.0e-9
    ./atmchange orbital_year=2019
    # use same GHG for land model
    ./xmlchange CCSM_CO2_PPMV=410.5
    #write out DAG
    ./atmchange atmosphere_dag_verbosity_level=5


    #specify land IC file
cat << EOF >> user_nl_elm
 finidat='/lustre/orion/cli115/world-shared/e3sm/inputdata/lnd/clm2/initdata/20220928.I2010CRUELM.ne1024pg2_ICOS10.elm.r.2016-08-01-00000.nc'

 hist_dov2xy = .true.,.true.
 hist_fincl2 = 'H2OSNO','SOILWATER_10CM','TG'
 hist_mfilt = 1,120
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','A'

EOF


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

    # Run type
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


    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.monthly_ne1024.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.50hourly_QcQiNcNi.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.50hourly_QrNrQmBm.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.6hourlyINST_ne30.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.6hourlyAVG_ne30.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.3hourlyAVG_ne120.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.3hourlyINST_ne120.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.3hourly_ne1024.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.hourly_2Dvars.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.ARM_sites_2D.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.ARM_sites_3D.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.monthly_cosp_ne1024.yaml" .
    cp ${SCREAMDOCS_ROOT}"/v1_output/scream_output.Cess.ACI_regions_2D.yaml" .

    ./atmchange output_yaml_files="./scream_output.Cess.monthly_ne1024.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.50hourly_QcQiNcNi.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.50hourly_QrNrQmBm.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.6hourlyINST_ne30.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.6hourlyAVG_ne30.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.3hourlyAVG_ne120.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.3hourlyINST_ne120.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.3hourly_ne1024.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.hourly_2Dvars.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.ARM_sites_2D.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.ARM_sites_3D.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.monthly_cosp_ne1024.yaml"
    ./atmchange output_yaml_files+="./scream_output.Cess.ACI_regions_2D.yaml"

    ./xmlchange --file env_run.xml --id SSTICE_DATA_FILENAME --val "/lustre/orion/cli115/world-shared/e3sm/inputdata/atm/cam/sst/sst_ostia_ukmo-l4_ghrsst_3600x7200_20190731_20200901_c20230913.nc"
    ./xmlchange --file env_run.xml --id  SSTICE_GRID_FILENAME --val "/lustre/orion/cli115/world-shared/e3sm/inputdata/ocn/docn7/domain.ocn.3600x7200.230522.nc"
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_ALIGN --val 2019
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_START --val 2019
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
    ./case.submit -a="-t ${WALLTIME} --mail-type=ALL --mail-user=terai1@llnl.gov"
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
