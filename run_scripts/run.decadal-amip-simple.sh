#!/bin/bash
set -e

branch="add-F20TR-SCREAMv1"
code_root=${HOME}/codes/scream/branches/${branch}
res=ne30pg2_EC30to60E2r2
compset=F20TR
machine=chrysalis
compiler=intel
pecount=5440x1
casename=${branch}.${res}.${compset}.${machine}_${compiler}.test_sstice1
caseroot=${HOME}/codes/scream/cases/${casename}
mkdir -p `dirname ${caseroot}`
${code_root}/cime/scripts/create_newcase \
    --case ${caseroot} \
    --res ${res} \
    --compset ${compset} \
    --pecount ${pecount} \
    --project e3sm \
    --walltime 01:00:00

cd ${caseroot}

# Extract input_data_dir for user edits to the namelist
input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

# Change run length
./xmlchange STOP_OPTION=nmonths,STOP_N=2
./xmlchange RESUBMIT=0
./xmlchange RUN_STARTDATE="1994-10-01"

# For big data
./xmlchange PIO_NETCDF_FORMAT="64bit_data"

# Namelist options for EAMxx
if [ "${res}" == "ne1024pg2_ne1024pg2" ]; then
    ./atmchange initial_conditions::Filename="\${DIN_LOC_ROOT}/atm/scream/init/screami_ne1024np4L128_era5-19941001-topoadjx6t_20240123.nc"
fi

# Turn on cosp and set frequency
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
# write out DAG
./atmchange atmosphere_dag_verbosity_level=5

if [ 1 -eq 0 ]; then  # BRHDEBUG
# Copy output stream yaml files
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
# Set output streams
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
fi


# Namelist options for ELM
cat << EOF > user_nl_elm
!fsurdat = $DIN_LOC_ROOT/lnd/clm2/surfdata_map/surfdata_ne30pg2_simyr1950_c210729.nc
!finidat = $DIN_LOC_ROOT/lnd/clm2/initdata/20210802.ICRUELM-1950.ne30pg2_EC30to60E2r2.elm.r.0051-01-01-00000.nc
EOF

    #specify land IC file
    if [ "${RESOLUTION}" == "ne1024pg2_ne1024pg2" ]; then
cat << EOF >> user_nl_elm
 ! Set input data paths
 finidat = "${DIN_LOC_ROOT}/lnd/clm2/initdata/20231226.I2010CRUELM.ne1024pg2_ICOS10.elm.r.1994-10-01-00000.nc"
 flanduse_timeseries = "${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_ne1024pg2_historical_simyr1990-2014_c240109.nc"
 fsurdat = "${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_ne1024pg2_historical_simyr1990-2014_c240109.nc"
EOF
    elif [ "${RESOLUTION}" == "ne30pg2_EC30to60E2r2" ]; then
cat << EOF >> user_nl_elm
 flanduse_timeseries = '${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_ne30np4.pg2_hist_simyr1850-2015_c210113.nc'
 fsurdat = '${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/surfdata_ne30pg2_simyr1850_c230417_with_TOP.nc'
EOF
    fi

cat << EOF >> user_nl_elm
 ! Override consistency checks since we know our surface data is inconsistent
 check_finidat_fsurdat_consistency = .false.
 check_finidat_year_consistency = .false.
 check_finidat_pct_consistency = .false.
 check_dynpft_consistency = .false.

 ! Make sure we do transient PFTs
 do_transient_pfts = .true.
 
 hist_dov2xy = .true.,.true.
 hist_fincl2 = 'H2OSNO','SOILWATER_10CM','TG'
 hist_mfilt = 1,120
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','A'

EOF

# Coupler settings; new surface flux scheme
cat << EOF >> user_nl_cpl
 ocn_surface_flux_scheme = 2
EOF


# Point to new SST forcing
if [ 1 -eq 0 ]; then
./xmlchange --file env_run.xml --id SSTICE_DATA_FILENAME --val "\${DIN_LOC_ROOT}/atm/cam/sst/sst_ostia_3600x7200_19941001_20151231_c20240108.nc"
./xmlchange --file env_run.xml --id SSTICE_GRID_FILENAME --val "\${DIN_LOC_ROOT}/ocn/docn7/domain.ocn.3600x7200.230522.nc"
./xmlchange --file env_run.xml --id SSTICE_YEAR_ALIGN --val 1994
./xmlchange --file env_run.xml --id SSTICE_YEAR_START --val 1994
./xmlchange --file env_run.xml --id SSTICE_YEAR_END --val 2015
fi

# Setup, build, run
./case.setup
./case.build
./case.submit
echo "${caseroot}"
