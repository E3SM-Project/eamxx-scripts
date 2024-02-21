# To use this test, first
#    git apply rpointer_test.patch

function run {
    echo "amb> CONTINUE_RUN is $contrun"
    echo "amb> crash $crash_ymd $crash_tod"
    echo "amb> stopn $stopn $unit"
    echo "amb> nresubmit $nresubmit"
    export amb_restart_crash_ymd=$crash_ymd
    export amb_restart_crash_tod=$crash_tod

    wcid=$account
    walltime="00:30:00"

    rm -rf $casename
    rm -f $rundir/../bld/cmake-bld/CMakeCache.txt

    ${repo}/cime/scripts/create_newcase \
             -case $casename -compset $compset -res $res \
             --machine $machine --compiler $compiler --project $wcid \
             --handle-preexisting-dirs u --walltime $walltime

    (cd $casename

     ./xmlchange GMAKE_J=$jmake
     ./xmlchange DEBUG=false

     ./xmlchange NTASKS=$(( $ncore * $nnode ))
     ./xmlchange NTHRDS=1
     ./xmlchange ROOTPE=0
     ./xmlchange MAX_MPITASKS_PER_NODE=$ncore
     ./xmlchange MAX_TASKS_PER_NODE=$ncore

     ./xmlchange STOP_OPTION=$unit
     ./xmlchange STOP_N=$stopn
     ./xmlchange REST_OPTION=$runit
     ./xmlchange REST_N=$restn
     ./xmlchange HIST_OPTION=$unit
     ./xmlchange HIST_N=$restn

     ./xmlchange JOB_QUEUE=$queue
     ./xmlchange SAVE_TIMING=FALSE

     ./xmlchange RESUBMIT=$nresubmit
     ./xmlchange CONTINUE_RUN=$contrun

     ./case.setup

     if [ $model == scream ]; then
         ./atmchange BfbHash=1 -b
     fi
     
     ./case.build
     ./case.submit
    )
}

function procphase {
    case $phase in
        1) args=(${phase1[*]}) ;;
        2) args=(${phase2[*]}) ;;
        3) args=(${phase3[*]}) ;;
    esac

    contrun=${args[0]}
    stopn=${args[1]}
    nresubmit=${args[2]}
    crash_ymd=${args[3]}
    crash_tod=${args[4]}
    echo $contrun $stopn $nresubmit $crash_ymd $crash_tod

    casename=rptr.$res.$compset.$compiler
    rundir=$outputbase/$casename/run
    echo rundir $rundir

    if [[ $phase == 1 ]]; then
        (cd $rundir; rm -f *log* *.nc *.bin rpointer.??? rpointer.???.prev)
    fi

    run

    if [[ $model == e3sm ]]; then
        grepline='^ nstep, te \s* $finstep '
        logfile='atm'
    else
        grepline='bfbhash> \s* $finstep '
        logfile='e3sm'
    fi

    (cd $rundir
     case $phase in
         1)
             while true; do
                 zgrep "rpointer>" e3sm.log.* > tmp.txt
                 grep crashing tmp.txt
                 if [ $? == 0 ]; then
                     break
                 else
                     echo "drive> waiting on 'crashing'"
                 fi
                 sleep 30
             done
             echo "drive> found 'crashing'"
             echo ">>> orig"; cat rpointer.???; echo ">>> prev"; cat rpointer.???.prev
             ;;
         2)
             while true; do
                 eval "zgrep \"${grepline}\" ${logfile}.log.*"
                 if [ $? == 0 ]; then
                     break
                 else
                     echo "drive> waiting on '$finstep'"
                 fi
                 sleep 30
             done
             echo "drive> found '$finstep'"
             ;;
         3)
             while true; do
                 eval "zgrep -l \"${grepline}\" ${logfile}.log.*" > tmp.txt
                 nmatch=(`wc -l tmp.txt`)
                 if [[ ${nmatch[0]} == 2 ]]; then
                     echo "drive> found '$finstep' twice:"
                     eval "zgrep \"${grepline}\" ${logfile}.log.*"
                     ls -ltrh ${logfile}.log.*
                     break
                 else
                     echo "drive> waiting on two matches for '$finstep'"
                 fi
                 sleep 30
             done
             ;;
     esac)
}

if [ $# -ne 2 ]; then
    echo "args: input-filename phase"
    exit
fi
input=$1
phase=$2

. $input

echo res $res
echo compset $compset
echo compiler $compiler
echo repo $repo
echo jmake $jmake
echo phase $phase

if [[ $phase == all ]]; then
    phase=1; procphase
    phase=2; procphase
    phase=3; procphase
else
    procphase
fi
