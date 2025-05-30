#!/bin/bash

# Set run options
#=============================================
resolution=ne1024pg2_r0125_oRRS18to6v3
compset=F2010-SCREAM-HR-DYAMOND2
checkout_date=2021015  #the date you *checked out* the code
branch=67ae74c          #actual git hash of branch to check out
run_descriptor=SCREAMv0 #will be SCREAMv0 for production run
repo=scream
machine=cori-knl
compiler=intel
stop_option="ndays"
stop_n="2"
rest_n="1"
walltime="32:00:00"
queue="regular"
debug_compile='FALSE'
date_string=`date`

# Setup processor layout
nnodes_atm=1536
nnodes_ocn=1536
nthreads=16
mpi_tasks_per_node=8
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
email_address="terai1@llnl.gov"

# Set flags specific to running this script
do_download=false
do_newcase=true
do_setup=true
do_build=true
do_submit=false

case_name=SCREAMv0.SCREAM-DY2.ne1024pg2.Feb11Case.20210115

# Set paths
#code_root=${HOME}/gitwork/scream/
code_root=${CSCRATCH}/E3SM_code/screamv020210115/
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
    ./xmlchange PIO_BUFFER_SIZE_LIMIT=134217728


    # Change the run type to branch
    ./xmlchange RUN_TYPE="branch"
    ./xmlchange RUN_REFDIR="/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127/run"
    ./xmlchange RUN_REFCASE="SCREAMv0.SCREAM-DY2.ne1024pg2.20201127"
    ./xmlchange RUN_REFDATE="2020-02-11"
    ./xmlchange GET_REFCASE="TRUE"

    
    # Edit CAM namelist to set dycore options for new grid
    cat <<EOF >> user_nl_eam

    !*** By default the model dumps hundreds of vars in h0. Don't do that. ***
    empty_htapes=.true.
    !*** Outputs for DYAMOND (note fincl can only go to 10) ***
    nhtfrq = 12,-6,-6,-6,-3,-3,-3,-3,-3,-3                  !output freq: 12 steps=15 mi, -6=6 hrs, -3=3hrs
    mfilt  = 96,4,4,4,8,8,8,8,8,8                            !new file freq: daily in all cases
    fincl1 = 'TVQ:I', 'TUQ:I', 'PS:I', 'PRECT:I'            !15 min 2d
    fincl2 = 'SHOC_TKE:A', 'WTHL_SEC:A',                    !bogenschutz priority 1,   6 hrly
    	     'WQW_SEC:A', 'OMEGAT:A', 'PS:A'
    fincl3 = 'OMEGAQ:A', 'SHOC_MIX:A', 'THL_SEC:A',         !bogenschutz priority 1-2 + Paul 6 hrly
    	     'W3:A', 'SNOWHLND:A'
    fincl4 = 'WTHV_SEC:A', 'TKH:A',                         !bogenschutz priority 2, 6 hrly
    	     'QW_SEC:A', 'QWTHL_SEC:A'
    fincl5 = 'P3_sed_CLDICE:I', 'P3_mtend_CLDICE:I',        !Hassan and Terai 3 hrly
    	     'RAINQM:I'
    fincl6 = 'DYN_PS:I', 'DYN_PNH:I'                        !Oksana 3 hrly
    fincl7 = 'VOR:I'                                        !Mark 3 hrly
    fincl8 = 'DIV:I'                                        !Mark and Hassan, 3 hrly
    fincl9 = 'DYN_U:I'                                      !Mark, priority 2 3 hrly
    fincl10 = 'DYN_V:I'                                     !Mark, priority 2 3 hrly
    !*** Rad Freq: 4x75sec=5 min which divides 15 min output freq ***
    iradsw = 4
    iradlw = 4

EOF

    cat <<EOF >> user_nl_elm


EOF
    
    # Finally, run setup
    ./case.setup

    # The run location is determined in the bowels of CIME
    # Symlink to that location from user-chosen $case_root (=current dir)
    ln -s `./xmlquery -value RUNDIR` run


    # This disables the logic that sets tprof_n and tprof_options internally.
    #./xmlchange --file env_run.xml TPROF_TOTAL=-1
    #echo "tprof_n = 1" >> user_nl_cpl
    #echo "tprof_option = 'nsteps'" >> user_nl_cpl
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

#==============================================
# IMPORTANT Note: CRT
./xmlchange GET_REFCASE="FALSE" # turn off GET_REFCASE so that the file doesn't try to copy over data again

# Here, you need to modify the rpointer files in the current run directory
# so that they have the correct dates that correspond to RUN_REFDATE

#==============================================


# Run
#=============================================
if [ "${do_submit}" == "true" ]; then
    cd ${case_root}
    ./case.submit --batch-args="--mail-type=ALL --mail-user=${email_address}"
fi

echo "Done working in ${case_root}"
