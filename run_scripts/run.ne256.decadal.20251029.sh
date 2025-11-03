#!/bin/bash -fe

# EAMxx template run script

main() {

do_fetch_code=true      
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

readonly MACHINE="pm-gpu"
readonly CHECKOUT="20251031"
readonly BRANCH="eamxx-v0.1.0"
readonly CHERRY=( )
readonly COMPILER="gnugpu"
readonly DEBUG_COMPILE=FALSE

readonly MPS_YES=FALSE
readonly MPS_NOS=8

readonly COMPSET="F20TR-SCREAMv1"
readonly RESOLUTION="ne256pg2_ne256pg2"

readonly CODE_ROOT="/pscratch/sd/t/terai/E3SM_code/eamxx-v0.1.0" #To change to path of own scratch
readonly PROJECT="e3sm"                                          #To change to own account

readonly TUNINGSET="spanc800.2xauto.acc150.n0032.test2.1"

readonly CASE_NAME=${RESOLUTION}.${COMPSET}.${CHECKOUT}.${TUNINGSET}

readonly CASE_ROOT="${SCRATCH}/e3sm_scratch/${MACHINE}/${CASE_NAME}"


readonly CASE_GROUP=""

readonly HIST_OPTION="never"
readonly HIST_N="1"

readonly MODEL_START_TYPE="initial"  # "initial", "continue", "branch", "hybrid"
readonly START_DATE="1994-10-01"     # "" for default, or explicit "0001-01-01"

readonly GET_REFCASE=false
readonly RUN_REFDIR=""
readonly RUN_REFCASE=""
readonly RUN_REFDATE=""   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
readonly CASE_RUN_DIR=${CASE_ROOT}/run
readonly NUM_PES=128

if [[ ${MPS_YES} == "TRUE" ]]; then
    readonly PEL_SIM=$(( MPS_NOS * NUM_PES ))
else
    readonly PEL_SIM=$(( 4 * NUM_PES ))
fi
readonly PELAYOUT=${PEL_SIM}"x1"

readonly WALLTIME="8:59:59"
readonly STOP_OPTION="nmonths"
readonly STOP_N="6"
readonly REST_OPTION="nmonths"
readonly REST_N="3"
readonly RESUBMIT="8"
readonly DO_SHORT_TERM_ARCHIVING=false
readonly OLD_EXECUTABLE=""

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

githash_eamxx=`git --git-dir ${CODE_ROOT}/.git rev-parse HEAD`

# Create case
create_newcase

# Copy script into case_script directory for provenance
copy_script

# Setup
case_setup

# Build
case_build

# Configure runtime options
runtime_options

# Submit
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
    local repo=E3SM

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

    if [[ ${MPS_YES} == "TRUE" ]]; then
        # MPS settings
        readonly MPS_NUMBER=${MPS_NOS}
    else
        # No MPS
        readonly MPS_NUMBER=4
    fi

    if [[ ${MPS_YES} == "TRUE" ]]; then
        ./xmlchange --file env_mach_pes.xml MAX_MPITASKS_PER_NODE="${MPS_NUMBER}"
        ./xmlquery --file env_mach_pes.xml MAX_MPITASKS_PER_NODE
    fi

    ./xmlchange --file env_mach_pes.xml NTHRDS="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ATM="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_LND="$(( 64 / MPS_NUMBER ))"
    ./xmlquery NTHRDS_LND
    ./xmlchange --file env_mach_pes.xml NTHRDS_ICE="$(( 64 / MPS_NUMBER ))"
    ./xmlquery NTHRDS_ICE
    ./xmlchange --file env_mach_pes.xml NTHRDS_OCN="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_ROF="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_CPL="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_GLC="1"
    ./xmlchange --file env_mach_pes.xml NTHRDS_WAV="1"

    ./xmlchange PIO_NETCDF_FORMAT="64bit_data"

    # Finally, run CIME case.setup
    ./case.setup --reset

    if [[ ${MPS_YES} == "TRUE" ]]; then
        # Use the mps wrapper
        ./xmlchange run_exe="${PSCRATCH}/mps-wrapper.sh `./xmlquery --value run_exe`"
        ./xmlquery run_exe
    fi

    # Save provenance invfo
    echo "branch hash for EAMxx: $githash_eamxx" > GIT_INFO.txt

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

    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # Set simulation start date
    if [ ! -z "${START_DATE}" ]; then
        ./xmlchange RUN_STARTDATE=${START_DATE}
    fi

    # Set temperature cut off in dycore threshold to 180K
    ./atmchange vtheta_thresh=180
    ./atmquery vtheta_thresh
    ./atmchange physics::atm_procs_list="mac_aero_mic,rrtmgp,cosp"
    ./atmchange physics::cosp::cosp_frequency_units="hours"
    ./atmchange physics::cosp::cosp_frequency=1

    # use GHG levels more appropriate for sim
    # Average from 19940101 - 20150101
    ./atmchange co2vmr=377.2e-6
    ./atmchange ch4vmr=1786.6e-9
    ./atmchange n2ovmr=318.6e-9
    ./atmchange orbital_year=-9999
    # use CO2 the same in land model
    ./xmlchange CCSM_CO2_PPMV=377.2
    #write out DAG
    #./atmchange atmosphere_dag_verbosity_level=5

    # new spa file, hey hoo
    ./atmchange spa_data_file="${input_data_dir}/atm/scream/init/spa_v3.LR.F2010.2011-2025.c_20240405.nc"
    ./atmchange spa_ccn_to_nc_factor=800.0

    # These are only needed if you want extra radiation diags (with aerosols, etc.)
    #./atmchange extra_clnclrsky_diag=true
    #./atmchange extra_clnsky_diag=true

# Specify land initial condition and surface datasets
cat << EOF >> user_nl_elm
  finidat = "${input_data_dir}/lnd/clm2/initdata/20240104.I2010CRUELM.ne256pg2.elm.r.1994-10-01-00000.nc"
  flanduse_timeseries = "${input_data_dir}/lnd/clm2/surfdata_map/landuse.timeseries_ne256pg2_hist_simyr1850-2015_c240131.nc"
  fsurdat = "${input_data_dir}/lnd/clm2/surfdata_map/surfdata_ne256pg2_simyr1850_c240131.nc"
EOF

# Additional settings for land
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
cat << EOF >> user_nl_mosart
  finidat_rtm = "${input_data_dir}/rof/mosart/initdata/20240104.I2010CRUELM.ne256pg2.mosart.r.1994-10-01-00000.nc"
EOF

# Coupler settings; new surface flux scheme
cat << EOF >> user_nl_cpl
  ocn_surface_flux_scheme = 2
EOF

# Turn off CICE completely IO completely (yuck)
cat << EOF >> user_nl_cice
 histfreq = 'x','x','x','x','x'
 histfreq_n = 0,0,0,0,0
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

cat <<EOF >> 1da.yaml
averaging_type: average
fields:
  physics_pg2:
    field_names:
    - snow_depth_land
max_snapshots_per_file: 10
filename_prefix: 1da
iotype: pnetcdf
output_control:
  frequency: 1
  frequency_units: ndays
restart:
  force_new_file: true
EOF

cat <<EOF >> 1hi.yaml
averaging_type: instant
fields:
  physics_pg2:
    field_names:
    - LW_flux_up_at_model_top
    - precip_total_surf_mass_flux
max_snapshots_per_file: 240
filename_prefix: 1hi
iotype: pnetcdf
output_control:
  frequency: 1
  frequency_units: nhours
restart:
  force_new_file: true
EOF

cat <<EOF >> 1ma_ne30pg2.yaml
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
    # COSP
    - isccp_ctptau
    - modis_ctptau
    - misr_cthtau
    - isccp_cldtot
max_snapshots_per_file: 1
filename_prefix: 1ma_ne30pg2
horiz_remap_file: \${DIN_LOC_ROOT}/atm/scream/maps/map_ne256pg2_to_ne30pg2_traave.20240206.nc
iotype: pnetcdf
output_control:
  frequency: 1
  frequency_units: nmonths
restart:
  force_new_file: true
EOF

cat <<EOF >> 3ha_ne30pg2.yaml
averaging_type: average
fields:
  physics_pg2:
    field_names:
    - T_mid
    - qv
    - qc
    - qr
    - qi
    - cldfrac_tot
    - cldfrac_liq
    - omega
    - U
    - V
    - pseudo_density
    - z_mid
    - surf_sens_flux
    - surf_evap
    - surface_upward_latent_heat_flux
    - ps
    - precip_liq_surf_mass_flux
    - precip_ice_surf_mass_flux
    - surf_mom_flux
    - surf_radiative_T
    - T_2m
    - sfc_flux_dir_nir
    - sfc_flux_dir_vis
    - sfc_flux_dif_nir
    - sfc_flux_dif_vis
    - sfc_flux_sw_net
    - sfc_flux_lw_dn
    - U_at_model_bot
    - V_at_model_bot
    - T_mid_at_model_bot
    - qv_at_model_bot
    - qc_at_model_bot
    - qi_at_model_bot
    - qr_at_model_bot
    - qm_at_model_bot
    - bm_at_model_bot
    - SW_flux_up_at_model_top
    - SW_flux_dn_at_model_top
    - LW_flux_up_at_model_top
    - SW_clrsky_flux_up_at_model_top
    - LW_clrsky_flux_up_at_model_top
    - SW_flux_dn_at_model_bot
    - SW_clrsky_flux_dn_at_model_bot
    - SW_flux_up_at_model_bot
    - SW_clrsky_flux_up_at_model_bot
    - LW_flux_dn_at_model_bot
    - LW_clrsky_flux_dn_at_model_bot
    - LW_flux_up_at_model_bot
    - LongwaveCloudForcing
    - ShortwaveCloudForcing
max_snapshots_per_file: 80
filename_prefix: 3ha_ne30pg2
horiz_remap_file: \${DIN_LOC_ROOT}/atm/scream/maps/map_ne256pg2_to_ne30pg2_traave.20240206.nc
iotype: pnetcdf
output_control:
  frequency: 3
  frequency_units: nhours
restart:
  force_new_file: true
EOF

cat <<EOF >> 1da_ne30pg2.yaml
averaging_type: average
fields:
  physics_pg2:
    field_names:
    - SeaLevelPressure
    - U_at_model_bot
    - V_at_model_bot
    - wind_speed_10m
    - U_at_925hPa
    - V_at_925hPa
    - U_at_850hPa
    - V_at_850hPa
    - U_at_700hPa
    - V_at_700hPa
    - U_at_500hPa
    - V_at_500hPa
    - U_at_300hPa
    - V_at_300hPa
    - ZonalVapFlux
    - MeridionalVapFlux
    - T_2m
    - qv_2m
    - T_mid_at_925hPa
    - T_mid_at_850hPa
    - T_mid_at_700hPa
    - T_mid_at_500hPa
    - T_mid_at_300hPa
    - RelativeHumidity_at_925hPa
    - RelativeHumidity_at_850hPa
    - RelativeHumidity_at_700hPa
    - RelativeHumidity_at_500hPa
    - RelativeHumidity_at_300hPa
    - z_mid_at_925hPa
    - z_mid_at_850hPa
    - z_mid_at_700hPa
    - z_mid_at_500hPa
    - z_mid_at_300hPa
    - SW_flux_up_at_model_top
    - SW_flux_dn_at_model_top
    - LW_flux_up_at_model_top
    - SW_clrsky_flux_up_at_model_top
    - SW_clrsky_flux_dn_at_model_top
    - LW_clrsky_flux_up_at_model_top
    - LiqWaterPath
    - precip_liq_surf_mass_flux
    - precip_ice_surf_mass_flux
    - IceWaterPath
    - cldfrac_tot_for_analysis
    - cldfrac_liq
    - omega_at_500hPa
    - omega_at_700hPa
    - omega_at_850hPa
    - PotentialTemperature_at_700hPa
    - PotentialTemperature_at_1000hPa
    - U_at_200hPa
    - V_at_200hPa
    - RelativeHumidity_at_200hPa
    - LiqNumberPath
    - IceNumberPath
    - eff_radius_qc
    - eff_radius_qi
    - surf_sens_flux
    - surf_evap
    - surface_upward_latent_heat_flux
    - ps
    - surf_mom_flux
    - surf_radiative_T
    - sfc_flux_dir_nir
    - sfc_flux_dir_vis
    - sfc_flux_dif_nir
    - sfc_flux_dif_vis
    - sfc_flux_sw_net
    - sfc_flux_lw_dn
    - T_mid_at_model_bot
    - qv_at_model_bot
    - qc_at_model_bot
    - qi_at_model_bot
    - qr_at_model_bot
    - qm_at_model_bot
    - bm_at_model_bot
    - SW_flux_dn_at_model_bot
    - SW_clrsky_flux_dn_at_model_bot
    - SW_flux_up_at_model_bot
    - SW_clrsky_flux_up_at_model_bot
    - LW_flux_dn_at_model_bot
    - LW_clrsky_flux_dn_at_model_bot
    - LW_flux_up_at_model_bot
    - LongwaveCloudForcing
    - ShortwaveCloudForcing
filename_prefix: 1da_ne30pg2
horiz_remap_file: \${DIN_LOC_ROOT}/atm/scream/maps/map_ne256pg2_to_ne30pg2_traave.20240206.nc
iotype: pnetcdf
max_snapshots_per_file: 10
output_control:
  frequency: 1
  frequency_units: ndays
restart:
  force_new_file: true
EOF

    ./atmchange output_yaml_files="./1da.yaml"
    ./atmchange output_yaml_files+="./1hi.yaml"
    ./atmchange output_yaml_files+="./1ma_ne30pg2.yaml"
    ./atmchange output_yaml_files+="./3ha_ne30pg2.yaml"
    ./atmchange output_yaml_files+="./1da_ne30pg2.yaml"

    # Point to new SST forcing (note 1994 start for F20TR)
    ./xmlchange --file env_run.xml --id SSTICE_DATA_FILENAME --val "${input_data_dir}/atm/cam/sst/sst_ostia_3600x7200_19940930_20151231_c20240125.nc"
    ./xmlchange --file env_run.xml --id SSTICE_GRID_FILENAME --val "${input_data_dir}/ocn/docn7/domain.ocn.3600x7200.230522.nc"
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_ALIGN --val 1994
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_START --val 1994
    ./xmlchange --file env_run.xml --id SSTICE_YEAR_END --val 2015


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
    # Add email settings and ignore some nodes via -x
    ./case.submit -a="-t ${WALLTIME} --mail-type=ALL --mail-user=terai1@llnl.gov"
    popd
}

#-----------------------------------------------------
copy_script() {

    echo $'\n----- Saving run script for provenance -----\n'
    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=$( basename -- "$0"; )
    local this_script_dir=$( dirname -- "$0"; )
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp "${this_script_dir}/${this_script_name}" ${script_provenance_dir}/${script_provenance_name}

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
