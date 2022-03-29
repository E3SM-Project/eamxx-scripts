#!/bin/bash

# Who to send email updates on run status to
email_address="bhillma@sandia.gov"

# Set run options
#resolution=ne512np4_360x720cru_oRRS15to5 #ne1024np4_360x720cru_oRRS15to5 #ne120_r0125_oRRS18to6v3 #ne4_ne4 #ne30_ne30 #ne1024np4_360x720cru_oRRS15to5
resolution=ne120np4_r0125_oRRS18to6v3 #ne1024np4_360x720cru_oRRS15to5 #ne120_r0125_oRRS18to6v3 #ne4_ne4 #ne30_ne30 #ne1024np4_360x720cru_oRRS15to5
compset=F2010-SCREAM-LR
#branch=68e0661e4 #  9bfb38267 # git hash to branch is most useful here
branch=aarondonahue/scream_p3_sa_input_from_f90 # git hash to branch is most useful here
repo=scream
machine=cori-knl
compiler=intel
stop_option="nsteps"
stop_n="1"
walltime="00:30:00"
queue="debug"

# Setup processor layout
nnodes=512 #384 #3072 #1536 #512 #128 #32 #1536 #12 #3072 #2048 #1536
nthreads=16
mpi_tasks_per_node=8
ntasks=$(expr ${nnodes} \* ${mpi_tasks_per_node})
total_tasks_per_node=$(expr ${mpi_tasks_per_node} \* ${nthreads})
pelayout=${nnodes}x${mpi_tasks_per_node}x${nthreads}

# Set flags specific to running this script
do_download=false
do_newcase=true
do_setup=true
do_build=true
do_submit=true

# Set paths
datestring=`date +"%Y%m%d"`
case_name=${branch}.${resolution}.${compset}.${machine}_${compiler}.${pelayout}.${datestring}
code_root=${HOME}/codes/${repo}/branches/${branch}
case_root=${SCRATCH}/${repo}/cases/${case_name}

# Download code
if [ "${do_download}" == "true" ]; then
    if [ ! -e ${code_root} ]; then
        echo "Cloning repository repo = $repo into branch = $branch under code_root = $code_root"
        cdir=`pwd`
        mkdir -p $code_root/

        # This will put repository, with all code, in directory $tag_name
        git clone git@github.com:E3SM-Project/${repo}.git $code_root
        
        # Setup git hooks
        rm -rf $code_root/.git/hooks
        git clone git@github.com:E3SM-Project/E3SM-Hooks.git $code_root/.git/hooks
        cd $code_root
        git config commit.template $code_root/.git/hooks/commit.template

        # Bring in all submodule components
        git submodule update --init

        # Check out desired branch
        git checkout ${branch}
        cd ${cdir}
    fi
fi

# Create case
if [ "${do_newcase}" == "true" ]; then
    ${code_root}/cime/scripts/create_newcase \
        --case ${case_root} \
        --compset ${compset} --res ${resolution} \
        --machine ${machine} --compiler ${compiler} \
        --queue ${queue} \
        --walltime ${walltime} || exit 1
fi

# Copy this script to case directory
cp -v `basename $0` ${case_root}/

# Setup
if [ "${do_setup}" == "true" ]; then
    cd ${case_root}

    # Set run length
    ./xmlchange STOP_OPTION=${stop_option},STOP_N=${stop_n}

    # Set processor layout
	./xmlchange NTASKS=${ntasks}
    ./xmlchange NTHRDS_ATM=${nthreads}
    ./xmlchange NTHRDS_LND=${nthreads}
	./xmlchange NTHRDS_CPL=1
    ./xmlchange MAX_MPITASKS_PER_NODE=${mpi_tasks_per_node}
    ./xmlchange MAX_TASKS_PER_NODE=${total_tasks_per_node}

    # Set PIO format, use PIO version 2, and increase PIO buffer size 
    ./xmlchange PIO_NETCDF_FORMAT="64bit_data"
    ./xmlchange PIO_BUFFER_SIZE_LIMIT=134217728

    if [ "${resolution}" == "ne120_r0125_oRRS18to6v3" ]; then #ne1024np4_360x720cru_oRRS15to5 #ne120_r0125_oRRS18to6v3 #ne4_ne4 #ne30_ne30 #ne1024np4_360x720cru_oRRS15to5
      cat <<EOF >> user_nl_eam
    ncdata = '/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme/cami_mam3_Linoz_0000-01-ne120np4_L72_c160318.nc'
EOF
    elif [ "${resolution}" == "ne512np4_360x720cru_oRRS15to5" ]; then
       cat <<EOF >> user_nl_eam
    ncdata = 'atm/cam/inic/homme/ifs_oper_T1279_2016080100_mod_subset_to_e3sm_ne512np4_topoadj_L128_c20200115.nc'
EOF
    fi

    # Edit CAM namelist to set dycore options for new grid
    cat <<EOF >> user_nl_eam

    ! Also write a new initial condition
    inithist = 'ENDOFRUN'

	! Always do radiation
	iradsw = 1
	iradlw = 1

    ! Outputs for initial conditions
    avgflag_pertape = 'A', 'I'
	nhtfrq = 0, 1
	mfilt = 1, 48
	fincl2 =
		"T_mid_inP3",
		"p_mid_inP3",
		"pseudo_density_inP3",
		"qc_inP3",
		"qv_inP3",
		"z_int_inP3",
		"inv_qc_relvar_inP3",
		"bm_inP3",
		"nc_inP3",
		"ni_inP3",
		"nr_inP3",
		"qi_inP3",
		"qm_inP3",
		"qr_inP3",
		"cldfrac_tot_inP3",
		"cldfrac_tot_inRAD",
		"T_prev_micro_step_inP3",
		"nc_activated_inP3",
		"nc_nuceat_tend_inP3",
		"ni_activated_inP3",
		"qv_prev_micro_step_inP3",
		"T_mid_inSHOC",
		"cldfrac_liq_inSHOC",
		"eddy_diff_mom_inSHOC",
		"horiz_winds_inSHOC",
		"host_dx_inSHOC",
		"host_dy_inSHOC",
		"omega_inSHOC",
		"p_int_inSHOC",
		"p_mid_inSHOC",
		"phis_inSHOC",
		"pseudo_density_inSHOC",
		"qc_inSHOC",
		"qv_inSHOC",
		"sgs_buoy_flux_inSHOC",
		"surf_latent_flux_inSHOC",
		"surf_sens_flux_inSHOC",
		"surf_mom_flux_inSHOC",
		"tke_inSHOC",
		"z_int_inSHOC",
		"z_mid_inSHOC",
		"T_mid_inRAD",
		"eff_radius_qc_inRAD",
		"eff_radius_qi_inRAD",
		"p_mid_inRAD",
		"p_int_inRAD",
		"qc_inRAD",
		"qi_inRAD",
		"surf_lw_flux_up_inRAD",
		"sfc_alb_dir_vis_inRAD",
		"sfc_alb_dif_vis_inRAD",
		"sfc_alb_dir_nir_inRAD",
		"sfc_alb_dif_nir_inRAD",
		"H2O_inRAD",
		"CO2_inRAD",
		"O3_inRAD",
		"N2O_inRAD",
		"CO_inRAD",
		"CH4_inRAD",
		"O2_inRAD",
		"N2_inRAD"
		"dp_inHOMME",
		"phi_int_inHOMME",
		"qv_inHOMME",
		"v_inHOMME",
		"vtheta_dp_inHOMME",
		"w_i_inHOMME",
		"ps_inHOMME"
EOF

    # Finally, run setup
    ./case.setup

    # Link to run directory
    ln -s `./xmlquery -value RUNDIR` run

    # Set file striping on run dir for writing large files
    ls -l run
	lfs setstripe -S 1m -c 32 run
    lfs getstripe run >& lfs_run.txt

    # This disables the logic that sets tprof_n and tprof_options internally.
    ./xmlchange --file env_run.xml TPROF_TOTAL=-1
    echo "tprof_n = 1" >> user_nl_cpl
    echo "tprof_option = 'nsteps'" >> user_nl_cpl

fi

# Build
if [ "${do_build}" == "true" ]; then
    cd ${case_root}
    ./case.build
fi

# Run
if [ "${do_submit}" == "true" ]; then
    cd ${case_root}
    ./case.submit --batch-args="--mail-type=ALL --mail-user=${email_address}"
fi
