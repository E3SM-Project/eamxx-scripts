#!/bin/bash

# Set run options
#=============================================
resolution=ne256pg2_r0125_oRRS18to6v3
compset=F2010-SCREAM-HR-DYAMOND2
checkout_date=20201112  #the date you *checked out* the code
branch=dyamond2         #actual name of branch to check out
run_descriptor=dyamond2 #will be SCREAMv0 for production run
repo=scream
machine=cori-knl
compiler=intel
stop_option="ndays"
stop_n="1"
rest_n="1"
walltime="2:00:00"
queue="regular"
debug_compile='FALSE'

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

case_name=${checkout_date}.${run_descriptor}.${compset}.${resolution}.${machine}.${pelayout}

# Set paths
code_root=${HOME}/gitwork/scream/
#code_root=${CSCRATCH}/E3SM_code/scream/
case_root=${CSCRATCH}/e3sm_scratch/cori-knl/${case_name}

######################################################################################
### USERS PROBABLY DON'T NEED TO CHANGE ANYTHING BELOW HERE EXCEPT user_nl_* FILES ###
######################################################################################

# Make directories created by this script world-readable:
#=============================================
umask 022

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
#=============================================
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
    nhtfrq = 3,3,3,3,-3,-3,-3,-3,-3,-3 !output freq: 3 steps=15 mi, -3=3hrs
    mfilt = 96,96,96,96,8,8,8,8,8,8 !new file freq: daily in all cases
    fincl1 = 'CLDLOW:I', 'CLDMED:I', 'CLDHGH:I', 'CLDTOT:I', 'TMCLDLIQ:I', 
             'TMCLDICE:I', 'TMRAINQM:I', 'TMCLDRIM:I', 'TMQ:I', 'CAPE:I', 'CIN:I'
    fincl2 = 'PS:I', 'TS:I', 'TREFHT:I', 'QREFHT:I', 'PRECT:I','PRECSL:I',
    	     'WINDSPD_10M:I', 'TAUX:I', 'TAUY:I', 'SHFLX:I', 'LHFLX:I'
    fincl3 = 'FSNTOA:I', 'FLNT:I','FLNTC:I','FSNTOAC:I', 'FSNS:I', 'FSDS:I', 
    	     'FLNS:I', 'FLDS:I'
    fincl4 = 'RH200:I',    'RH500:I',    'RH700:I',    'RH850:I',
	     'OMEGA200:I', 'OMEGA500:I', 'OMEGA700:I', 'OMEGA850:I', 
             'Z200:I',     'Z500:I',     'Z700:I',     'Z850:I'
    !*** 3 hrly (mostly 3d) variables below here ***
    fincl5 = 'PS:I', 'PSL:I', 'TMNUMLIQ:I', 'TMNUMICE:I', 'TMNUMRAI:I'
    fincl6 = 'U:I', 'V:I'
    fincl7 = 'T:I', 'Q:I', 
    fincl8 = 'CLDLIQ:I', 'CLDICE:I'
    fincl9 = 'CLOUD:I','OMEGA:I'
    fincl10= 'EMIS:I', 'TOT_ICLD_VISTAU:I'
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
    ls -l run > /dev/null #just using this command as a check that run dir exists
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
