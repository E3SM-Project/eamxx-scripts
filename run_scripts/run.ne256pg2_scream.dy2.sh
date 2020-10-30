#!/bin/bash

# Set run options
#=============================================
resolution=ne256pg2_r0125_oRRS18to6v3
compset=F2010-SCREAM-HR-DYAMOND2
branch=autoconv_25perc
repo=scream
machine=cori-knl
compiler=intel
stop_option="ndays"
stop_n="1"
rest_n="1"
walltime="2:00:00"
queue="regular"

# Setup processor layout
nnodes_atm=1024
nnodes_ocn=1024
nthreads=8
mpi_tasks_per_node=16
ntasks_atm=$(expr ${nnodes_atm} \* ${mpi_tasks_per_node})
ntasks_ocn=$(expr ${nnodes_ocn} \* ${mpi_tasks_per_node})
total_tasks_per_node=$(expr ${mpi_tasks_per_node} \* ${nthreads})
if [ ${ntasks_ocn} -ne ${ntasks_atm} ]; then
    nnodes=$(expr ${nnodes_atm} + ${nnodes_ocn})
else
    nnodes=${nnodes_atm}
fi
pelayout=${nnodes}x${mpi_tasks_per_node}x${nthreads}

# Who to send email updates on run status to
email_address="caldwell19@llnl.gov"

# Set flags specific to running this script
do_download=false
do_newcase=true
do_setup=true
do_build=true
do_submit=true

debug_compile='FALSE'

# Set paths
datestring=`date +"%Y%m%d-%H"`
case_name=${branch}.${resolution}.${compset}.${machine}_${compiler}.${pelayout}.DY2_Oct6.${datestring}
code_root=${HOME}/gitwork/scream/
case_root=${CSCRATCH}/E3SM_runs/${case_name}

######################################################################################
### USERS PROBABLY DON'T NEED TO CHANGE ANYTHING BELOW HERE EXCEPT user_nl_* FILES ###
######################################################################################

# Download code
#=============================================
if [ "${do_download}" == "true" ]; then

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
    git submodule update --init --recursive

    # Check out desired branch
    git checkout ${branch}
    cd ${cdir}

fi

# Create case
#=============================================
if [ "${do_newcase}" == "true" ]; then
    ${code_root}/cime/scripts/create_newcase \
        --case ${case_root} \
        --compset ${compset} --res ${resolution} \
        --machine ${machine} --compiler ${compiler} \
        --queue ${queue} \
        --walltime ${walltime} 
fi

# Copy this script to case directory
cp -v `basename $0` ${case_root}/

# Setup
#=============================================
if [ "${do_setup}" == "true" ]; then
    cd ${case_root}

    # Set run length
    ./xmlchange STOP_OPTION=${stop_option},STOP_N=${stop_n}
    ./xmlchange REST_N=${rest_n}

    # Set processor layout
    if [ ${ntasks_ocn} -ne ${ntasks_atm} ]; then
        ./xmlchange NTASKS=${ntasks_ocn}
        ./xmlchange NTASKS_ATM=${ntasks_atm}
        ./xmlchange ROOTPE_ATM=${ntasks_ocn}
    else
        ./xmlchange NTASKS=${ntasks_atm}
    fi
    ./xmlchange NTHRDS_ATM=${nthreads}
    ./xmlchange MAX_MPITASKS_PER_NODE=${mpi_tasks_per_node}
    ./xmlchange MAX_TASKS_PER_NODE=${total_tasks_per_node}

    # Flag for debug compile
    ./xmlchange --id DEBUG --val ${debug_compile}
    
    # Set PIO format, use PIO version 2, and increase PIO buffer size 
    ./xmlchange PIO_NETCDF_FORMAT="64bit_data"
    ./xmlchange PIO_VERSION="2"
    ./xmlchange PIO_BUFFER_SIZE_LIMIT=64200000
    
    # Edit CAM namelist to set dycore options for new grid
    cat <<EOF >> user_nl_eam

    !*** By default the model dumps hundreds of vars in h0. Don't do that. ***
    empty_htapes=.true.
    !*** Outputs for DYAMOND (note fincl can only go to 10) ***
    nhtfrq = 3,3,3,3,3,-3,-3,-3,-3,-3 !output freq: 3 steps=15 mi, -3=3hrs
    mfilt = 96,96,96,96,96,8,8,8,8,8 !new file freq: daily in all cases
    fincl1 = 'CLDLOW', 'CLDMED', 'CLDHGH', 'CLDTOT', 
             'TMCLDLIQ', 'TMCLDICE', 'TMRAINQM', 'TMCLDRIM', 'TMQ' 
    fincl2 = 'PS', 'TS', 'TREFHT', 'QREFHT',  
             'SHFLX', 'LHFLX', 
	     'PRECT','PRECSL', 'QFLX'  
             !'CAPE', 'CIN', 'V10'
    fincl3 = 'FSNTOA', 'FLNT','FLNTC','FSNTOAC', 'FSNS', 'FSDS', 'FLNS', 'FLDS'
    fincl4 = 'U10', 'TAUX', 'TAUY' !note U10 is wind speed. We plan to change this.
    fincl5 = 'T200',     'T500',     'T700',     'T850', !remove in favor of RH?
             'Q200',     'Q500',     'Q700',     'Q850', !remove in favor of RH?
	     'OMEGA200', 'OMEGA500', 'OMEGA700', 'OMEGA850', 
             'Z200',     'Z500',     'Z700',     'Z850'
    !*** 3d variables below here ***
    fincl6 = 'PS:I'
    fincl7 = 'U:I', 'V:I', 'OMEGA:I'
    fincl8 = 'T:I', 'Q:I', 
    fincl9 = 'CLDLIQ:I', 'CLDICE:I'
    fincl10 = 'NUMICE:I','NUMLIQ:I' !remove in favor of something else?
    !*** Radiation must be called every output timestep ***
    iradsw = 3
    iradlw = 3

EOF

    cat <<EOF >> user_nl_elm


EOF
    
    # Finally, run setup
    ./case.setup

    # The run location is determined in the bowels of CIME
    # Symlink to that location from user-chosen $case_root (=current dir)
    ln -s `./xmlquery -value RUNDIR` run


    # This disables the logic that sets tprof_n and tprof_options internally.
    ./xmlchange --file env_run.xml TPROF_TOTAL=-1
    echo "tprof_n = 1" >> user_nl_cpl
    echo "tprof_option = 'nsteps'" >> user_nl_cpl
fi

# Build
#=============================================
if [ "${do_build}" == "true" ]; then
    cd ${case_root}
    ./case.build

    # Set file striping on run dir for writing large files
    ls -l run
    lfs setstripe -S 1m -c 64 run
    lfs getstripe run >& lfs_run.txt
fi

# Run
#=============================================
if [ "${do_submit}" == "true" ]; then
    cd ${case_root}
    ./case.submit --batch-args="--mail-type=ALL --mail-user=${email_address}"
fi

echo "Done working in ${case_root}"
