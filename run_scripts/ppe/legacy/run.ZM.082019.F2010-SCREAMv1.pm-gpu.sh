#!/bin/bash -fe

# EAMxx PPE run script: ZM + SPC, L128v4, Aug-2019 start.
# One file for every resolution (ne32, ne128, ne256) and both modes (control, ensemble member).
#
# Control (builds and submits):
#   ./run.ZM.082019.F2010-SCREAMv1.pm-gpu.sh ne128
#   do_case_build=false ./run.ZM.082019.F2010-SCREAMv1.pm-gpu.sh ne32     # reuse the control's own build
#   Q=debug WALLTIME=00:30:00 ./run.ZM.082019.F2010-SCREAMv1.pm-gpu.sh ne32  # quick test
#
# Ensemble members (create_ensembles.py sets MEMBER_ID, the 19 PPE vars, do_case_build=false, OLD_EXECUTABLE):
#   RES=ne128 python create_ensembles.py --template ./run.ZM.082019.F2010-SCREAMv1.pm-gpu.sh \
#       --json ./append_nested_ne128_256_lhs.json --old-exe <ne128 control case>/build --start 0 --count 32

# create_ensembles.py runs "bash <template>", which ignores the shebang flags, so set them explicitly.
set -fe

# PPE parameters: env var set by create_ensembles.py : namelist key. Must match PARAMS in create_ensembles.py.
readonly PPE_PARAMS=(
	SHOC_THL2TUNE:shoc::thl2tune
	SHOC_QW2TUNE:shoc::qw2tune
	SHOC_LENGTH_FAC:shoc::length_fac
	SHOC_C_DIAG_3RD_MOM:shoc::c_diag_3rd_mom
	SHOC_COEFF_KH:shoc::coeff_kh
	SHOC_COEFF_KM:shoc::coeff_km
	SHOC_LAMBDA_LOW:shoc::lambda_low
	SHOC_LAMBDA_HIGH:shoc::lambda_high
	P3_SPA_CCN_TO_NC_FACTOR:p3::spa_ccn_to_nc_factor
	P3_CLDLIQ_TO_ICE_COLLECTION_FACTOR:p3::cldliq_to_ice_collection_factor
	P3_RAIN_TO_ICE_COLLECTION_FACTOR:p3::rain_to_ice_collection_factor
	P3_ACCRETION_PREFACTOR:p3::accretion_prefactor
	P3_DEPOSITION_NUCLEATION_EXPONENT:p3::deposition_nucleation_exponent
	P3_MAX_TOTAL_NI:p3::max_total_ni
	P3_ICE_SEDIMENTATION_FACTOR:p3::ice_sedimentation_factor
	P3_RAIN_SELFCOLLECTION_BREAKUP_DIAMETER:p3::rain_selfcollection_breakup_diameter
	P3_AUTOCONVERSION_PREFACTOR:p3::autoconversion_prefactor
	P3_AUTOCONVERSION_QC_EXPONENT:p3::autoconversion_qc_exponent
	P3_AUTOCONVERSION_RADIUS:p3::autoconversion_radius
)

main() {

	# Resolution: first argument, or RES from the environment (create_ensembles.py passes no arguments)
	RES="${1:-${RES:-}}"
	# create_ensembles.py always sets MEMBER_ID, which selects ensemble mode
	MEMBER_ID="${MEMBER_ID:-}"
	if [ -n "${MEMBER_ID}" ]; then readonly MODE=ensemble; else readonly MODE=control; fi

	# Control builds and submits; members reuse the control build and are submitted later as a bundle
	if [ "${MODE}" == "control" ]; then
		do_case_build="${do_case_build:-true}"
		do_case_submit="${do_case_submit:-true}"
	else
		do_case_build="${do_case_build:-false}"
		do_case_submit="${do_case_submit:-false}"
	fi
	do_create_newcase="${do_create_newcase:-true}"
	do_case_setup="${do_case_setup:-true}"

	readonly MACHINE="pm-gpu"
	readonly CHECKOUT="d5090c849a40dbe727c2f240e917b5c64de75f97"
	readonly COMPILER="gnugpu"
	readonly DEBUG_COMPILE=FALSE

	# Simulation
	readonly COMPSET="F2010-SCREAMv1"
	resolution_settings

	readonly CODE_ROOT="/pscratch/sd/b/beydoun/e3sm_repo_10132026/E3SM"
	readonly PROJECT="e3sm"
	githash_eamxx=`git --git-dir ${CODE_ROOT}/.git rev-parse HEAD`

	# Members are named after their control case and live in a sibling *_ens directory
	readonly NODES=$(( ${PELAYOUT%x*} / 4 ))
	readonly CONTROL_NAME=PPEensemble_${NODES}node.${RESOLUTION}.${COMPSET}.${CHECKOUT}
	readonly CONTROL_ROOT_BASE="${SCRATCH}/e3sm_scratch/pm-gpu/${RES}_ppe_v4_config_new_workflow"
	if [ "${MODE}" == "control" ]; then
		readonly CASE_NAME=${CONTROL_NAME}
		readonly CASE_ROOT="${CONTROL_ROOT_BASE}/${CASE_NAME}"
	else
		readonly CASE_NAME=${CONTROL_NAME}.${MEMBER_ID}
		readonly CASE_ROOT="${CONTROL_ROOT_BASE}_ens/${CASE_NAME}"
	fi
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
	readonly Q="${Q:-regular}"
	readonly WALLTIME="${WALLTIME:-${DEFAULT_WALLTIME}}"
	readonly STOP_OPTION="nmonths"
	readonly STOP_N="13"
	readonly REST_OPTION="nmonths"
	readonly REST_N="3"
	readonly RESUBMIT="0"
	readonly DO_SHORT_TERM_ARCHIVING=false
	# Build directory to reuse (members: the control's build dir, set by create_ensembles.py --old-exe)
	OLD_EXECUTABLE="${OLD_EXECUTABLE:-}"

	echo "RES=${RES} MODE=${MODE} CASE_NAME=${CASE_NAME} OLD_EXECUTABLE=${OLD_EXECUTABLE}"

	# --- Now, do the work ---
	# Make directories created by this script world-readable
	umask 022
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
	case_submit
	# All done
	echo $'\n----- All done -----\n'
}

#-----------------------------------------------------
# The only place the resolutions differ. Paths without a leading / are relative to DIN_LOC_ROOT.
resolution_settings() {
	local out_maps=/global/cfs/cdirs/e3sm/beydoun/autotune_utils
	case "${RES}" in
	ne32)
		RESOLUTION=ne32pg2_ne32pg2
		PELAYOUT=4x1
		DEFAULT_WALLTIME=06:00:00
		ATM_NCPL=""
		SPC_REMAP_FILE=atm/scream/maps/map_ne30pg2_to_ne32pg2_traave.20251124.nc
		OUT_REMAP_FILE=""   # output stays on native ne32pg2
		LND_FINIDAT=""
		;;
	ne128)
		RESOLUTION=ne128pg2_ne128pg2
		PELAYOUT=64x1
		DEFAULT_WALLTIME=30:00:00
		ATM_NCPL=96         # dtime=15min for ne128
		# no spc_remap_file default for ne128 at this checkout; same map SPA uses for ne128
		SPC_REMAP_FILE=atm/scream/maps/map_ne30pg2_to_ne128pg2_traave.20251124.nc
		OUT_REMAP_FILE=${out_maps}/map_ne128pg2_to_ne32pg2.nc
		LND_FINIDAT=lnd/clm2/initdata_map/ne128pg2.elm.r.2013-08-01-00000.64bit.nc
		;;
	ne256)
		RESOLUTION=ne256pg2_ne256pg2
		PELAYOUT=256x1
		DEFAULT_WALLTIME=30:00:00
		ATM_NCPL=""
		SPC_REMAP_FILE=atm/scream/maps/map_ne30pg2_to_ne256pg2_20231201.nc
		OUT_REMAP_FILE=${out_maps}/map_ne256pg2_to_ne32pg2.nc
		LND_FINIDAT=lnd/clm2/initdata/20240104.I2010CRUELM.ne256pg2.elm.r.2016-08-01-00000.nc
		;;
	*)
		echo "ERROR: resolution must be ne32, ne128 or ne256 (got '${RES}')."
		echo "       Pass it as the first argument, or as RES=... in front of create_ensembles.py."
		exit 1
		;;
	esac
	# Aug 1 initial condition on the L128v4 grid
	IC_FILE=/global/cfs/projectdirs/e3sm/whannah/HICCUP/HICCUP.atm_era5.2019-08-01.${RES}np4.L128v4.20261006.nc
}

######################################################
### Most users won't need to change anything below ###
######################################################
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
    # Optional arguments
    if [ ! -z "${PROJECT}" ]; then
      args="${args} --project ${PROJECT}"
    fi
    if [ ! -z "${CASE_GROUP}" ]; then
      args="${args} --case-group ${CASE_GROUP}"
    fi
    if [ ! -z "${QUEUE}" ]; then
      args="${args} --queue ${QUEUE}"
    fi
	# run create_newcase
    if ! ${CODE_ROOT}/cime/scripts/create_newcase ${args}; then
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
    # Extract input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`
	# switch to new vertical grid (also comes with new IC file)
	./xmlchange EAMXX_VGRID=L128v4
	# Set threading
    ./xmlchange --file env_mach_pes.xml NTHRDS="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ATM="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_LND="16"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ICE="16"
    ./xmlchange --file env_mach_pes.xml NTHRDS_OCN="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ROF="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_CPL="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_GLC="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_WAV="1"
    ./xmlchange PIO_NETCDF_FORMAT="64bit_data"
    # Finally, run CIME case.setup
    ./case.setup --reset
    # Save provenance info
    echo "branch hash for EAMxx: $githash_eamxx" > GIT_INFO.txt
    popd
}
#-----------------------------------------------------
case_build() {
    pushd ${CASE_SCRIPTS_DIR}
    # do_case_build = false
    if [ "${do_case_build,,}" != "true" ]; then
		echo $'\n----- case_build -----\n'
		if [ -n "${OLD_EXECUTABLE}" ]; then
		    # Reuse another case's build directory (members: the control build of the same resolution)
		    if [ ! -x "${OLD_EXECUTABLE}/e3sm.exe" ]; then
			echo "ERROR: OLD_EXECUTABLE must be a build dir containing e3sm.exe; got '${OLD_EXECUTABLE}'"
			exit 297
		    fi
		    if [[ "${OLD_EXECUTABLE}" != *"${RESOLUTION}"* ]]; then
			echo "ERROR: OLD_EXECUTABLE '${OLD_EXECUTABLE}' does not look like a ${RESOLUTION} build (RES=${RES})"
			exit 297
		    fi
		    ./xmlchange EXEROOT=${OLD_EXECUTABLE}
		    ./case.setup --reset
		elif [ -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
		    echo 'Skipping build because $do_case_build = '${do_case_build}
		else
		    echo 'ERROR: $do_case_build = '${do_case_build}' but no executable exists for this case.'
		    exit 297
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
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`
    # Set simulation start date
    if [ ! -z "${START_DATE}" ]; then
	./xmlchange RUN_STARTDATE=${START_DATE}
    fi
	# override the default IC file to be consistent with the Aug 1 start date
	./atmchange initial_conditions::filename="${IC_FILE}"

    ./atmchange -b vertical_coordinate_filename="${input_data_dir}/atm/scream/init/vertical_coordinates_L128v4_c20260820.nc"
	# Add ZM and SPC
	./atmchange -b physics::atm_procs_list="zm,mac_aero_mic,spc,rrtmgp"
    # use updated spa file (ne30pg2, remapped to the model grid)
    ./atmchange -b spa_data_file="${input_data_dir}/atm/scream/init/spa_from_v3LRamip_CMIP6_2010_clim_ne30pg2_c20260928.nc"
    ./atmchange -b spc_remap_file="${input_data_dir}/${SPC_REMAP_FILE}"
	# various other parameter updates
    ./atmchange -b set_cld_frac_r_to_one=True
	# disable split-form - despite being helpful at high-resolution, it can be problematic for low-res and we don't understand why
	# ./atmchange -b pgrad_correction=0
	# ./atmchange -b theta_advect_form=2
    if [ -n "${ATM_NCPL}" ]; then
	./xmlchange ATM_NCPL=${ATM_NCPL}
    fi
    #set sst inputs
    ./xmlchange --file env_run.xml --id SSTICE_DATA_FILENAME --val "${input_data_dir}/atm/cam/sst/sst_ostia_ukmo-l4_ghrsst_3600x7200_20190731_20210309_c20240506.nc"
    ./xmlchange --file env_run.xml --id SSTICE_GRID_FILENAME --val "${input_data_dir}/ocn/docn7/domain.ocn.3600x7200.230522.nc"
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_ALIGN --val 2019
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_START --val 2019
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_END --val 2021
    # PPE parameters (ensemble members only)
    ppe_parameters
    #./atmchange BfbHash=1
    #./atmchange --all internal_diagnostics_level=1 atmosphere_processes::internal_diagnostics_level=1
    #./atmchange ANY::internal_diagnostics_level=1
    # land history and IC (guarded so re-running a case doesn't duplicate entries)
	if ! grep -q hist_fincl2 user_nl_elm; then
cat << EOF >> user_nl_elm
 hist_dov2xy = .true.,.true.
 hist_fincl2 = 'H2OSNO','SOILWATER_10CM','TG'
 hist_mfilt = 1,120
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','A'
EOF
	    if [ -n "${LND_FINIDAT}" ]; then
		echo " finidat=\"${input_data_dir}/${LND_FINIDAT}\"" >> user_nl_elm
	    fi
	fi
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

    # Output: 1ma and 3ha are remapped to ne32pg2 above ne32 (hence the _ne32pg2 names); 51hi stays native.
    # Files are written with ">" so re-running a case doesn't duplicate YAML keys.
    local remap=""
    if [ -n "${OUT_REMAP_FILE}" ]; then
	remap=$'\n'"horiz_remap_file: ${OUT_REMAP_FILE}"
    fi
cat <<EOF > 1ma_ne32pg2.yaml
averaging_type: average
fields:
  physics_pg2:
    field_names:
    # 3D fields
    - T_mid
    - qv
    - RelativeHumidity
    - qc
    - qi
    - qr
    - qm
    - nc
    - ni
    - nr
    - cldfrac_tot_for_analysis
    - cldfrac_ice_for_analysis
    - cldfrac_liq
    - omega
    - U
    - V
    - z_mid
    - p_mid
    - tke
    # 2D fields
    - SW_flux_up_at_model_top
    - SW_flux_dn_at_model_top
    - LW_flux_up_at_model_top
    - SW_clrsky_flux_up_at_model_top
    - SW_clrsky_flux_dn_at_model_top
    - LW_clrsky_flux_up_at_model_top
    - SW_flux_up_at_model_bot
    - SW_flux_dn_at_model_bot
    - LW_flux_up_at_model_bot
    - LW_flux_dn_at_model_bot
    - SW_clrsky_flux_up_at_model_bot
    - SW_clrsky_flux_dn_at_model_bot
    - LW_clrsky_flux_dn_at_model_bot
    - ShortwaveCloudForcing
    - LongwaveCloudForcing
    - ps
    - SeaLevelPressure
    - T_2m
    - qv_2m
    - surf_radiative_T
    - VapWaterPath
    - IceWaterPath
    - LiqWaterPath
    - RainWaterPath
    - ZonalVapFlux
    - MeridionalVapFlux
    - surf_evap
    - surf_sens_flux
    - surface_upward_latent_heat_flux
    - precip_liq_surf_mass_flux
    - precip_ice_surf_mass_flux
    - landfrac
    - ocnfrac
    - PotentialTemperature_at_700hPa
    - PotentialTemperature_at_850hPa
    - PotentialTemperature_at_1000hPa
    - PotentialTemperature_at_2m_above_surface
    - omega_at_500hPa
    - omega_at_700hPa
    - omega_at_850hPa
    - RelativeHumidity_at_700hPa
    - RelativeHumidity_at_1000hPa
    - RelativeHumidity_at_2m_above_surface
    - wind_speed_10m
    - z_mid_at_700hPa
    - z_mid_at_1000hPa
    - T_mid_at_850hPa
    - T_mid_at_700hPa
    # For SST advection
    - U_at_10m_above_surface
    - V_at_10m_above_surface
max_snapshots_per_file: 1
filename_prefix: 1ma_ne32pg2${remap}
iotype: pnetcdf
output_control:
  frequency: 1
  frequency_units: nmonths
restart:
  force_new_file: true
EOF
cat <<EOF > 3ha_ne32pg2.yaml
averaging_type: average
fields:
  physics_pg2:
    field_names:
    - precip_total_surf_mass_flux
    - U_at_850hPa
    - V_at_850hPa
    - LW_flux_up_at_model_top
max_snapshots_per_file: 8
filename_prefix: 3ha_ne32pg2${remap}
iotype: pnetcdf
output_control:
  frequency: 3
  frequency_units: nhours
restart:
  force_new_file: false
EOF
cat <<EOF > 51hi.yaml
averaging_type: instant
fields:
  physics_pg2:
    field_names:
    - LW_flux_up_at_model_top
    - SW_flux_up_at_model_top
    - precip_total_surf_mass_flux
    - T_2m
    - VapWaterPath
    - IceWaterPath
    - LiqWaterPath
    - RainWaterPath
max_snapshots_per_file: 8
filename_prefix: 1hi
iotype: pnetcdf
output_control:
  frequency: 51
  frequency_units: nhours
restart:
  force_new_file: false
EOF
    ./atmchange output_yaml_files="./1ma_ne32pg2.yaml"
    ./atmchange output_yaml_files+="./3ha_ne32pg2.yaml"
    ./atmchange output_yaml_files+="./51hi.yaml"
    # Log the process list and the 19 PPE values actually in the namelist (members: ensemble_logs/<member>.log)
    local p keys=()
    for p in "${PPE_PARAMS[@]}"; do keys+=("${p#*:}"); done
    ./atmquery physics::atm_procs_list "${keys[@]}"
    popd
}
#-----------------------------------------------------
ppe_parameters() {
    if [ "${MODE}" != "ensemble" ]; then
	return
    fi
    local p env changes=()
    for p in "${PPE_PARAMS[@]}"; do
	env=${p%%:*}
	if [ -z "${!env:-}" ]; then
	    echo "ERROR: Missing ${env}"
	    exit 1
	fi
	changes+=("${p#*:}=${!env}")
    done
    # No -b: atmchange validates every name now, so a bad one fails member creation, not bundle submit
    ./atmchange "${changes[@]}"
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
    ./case.submit -a="-t ${WALLTIME} --qos=${Q} --mail-type=ALL --mail-user=beydoun1@llnl.gov"
    popd
}
#-----------------------------------------------------
copy_script() {
    echo $'\n----- Saving run script for provenance -----\n'
    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=`basename $0`
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    # copy from $0 (absolute path when run by create_ensembles.py), not the cwd
    cp -vp "$0" ${script_provenance_dir}/${script_provenance_name}
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
main "$@"
