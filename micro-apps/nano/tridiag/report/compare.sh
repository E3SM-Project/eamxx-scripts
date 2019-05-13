# Compare scream tridiag performance test against scream-docs standalone tridiag
# nano-app. This is to assure that the scream version is getting the right build
# flags and gives the right performance.

machine=$1 # v100, skx, knl
nrep=11

standalone=./standalone
scream=./tridiag

prefix=""
if [ $machine == "knl" ]; then
    export OMP_NUM_THREADS=136 # or 272
    export OMP_PROC_BIND=false
    export KMP_AFFINITY=balanced,granularity=fine
    prefix=" numactl -i 1"
elif [ $machine == "v100" ]; then
    prefix=""
elif [ $machine == "ws" ]; then
    export OMP_NUM_THREADS=4
    export OMP_PROC_BIND=false
    export KMP_AFFINITY=balanced,granularity=fine
    prefix=""    
else
    export OMP_NUM_THREADS=96
    export OMP_PROC_BIND=false
    export KMP_AFFINITY=balanced,granularity=fine
    prefix=""
fi

if [ $machine == "v100" ]; then
    for trial in $(seq $nrep); do
        # Compare against cusparse on cusparse optimal problem
        echo ">>> fast"
        for nlev in 72 128; do
            for nprob in 96 136 192 272 384 768 1632 2176 4096 8192; do
                for nw in 1 2 4 8; do
                    $standalone-$machine -np $nprob -nr $nlev -nw $nw -nc 1 -m cr_a1x1
                    $standalone-$machine -np $nprob -nr $nlev -nw $nw -nc 1 -m cr_a1x1p
                    $standalone-$machine -np $nprob -nr $nlev -nw $nw -nc 1 -m cr_a1xm
                    $standalone-$machine -np $nprob -nr $nlev -nw $nw -nc 1 -m cr_amxm
                    $standalone-$machine -np $nprob -nr $nlev -nw $nw -nc 1 -m thomas
                    $scream-$machine -np $nprob -nr $nlev -nc 1 -nw $nw -m cr
                    $scream-$machine -np $nprob -nr $nlev -nc 1 -nw $nw -m thomas
                done
            done
        done
        # Multiple RHS
        echo ">>> nrhs"
        for nrhs in 1 2 13 16 43; do
            for nlev in 72 128; do
                for nprob in 96 136 192 272 384 768 1632 2176; do
                    for nw in 1 2 4 8 16; do
                        $standalone-$machine -np $nprob -nr $nlev -nw $nw -nc $nrhs -m cr_a1xm
                        $standalone-$machine -np $nprob -nr $nlev -nw $nw -nc $nrhs -m thomas
                        $scream-$machine -np $nprob -nr $nlev -nc $nrhs -nw $nw -m cr -1a
                        $scream-$machine -np $nprob -nr $nlev -nc $nrhs -nw $nw -m thomas -1a
                    done
                done
            done
        done
        # Homme problem
        echo ">>> homme"
        for nrhs in 4 8 10 16 20; do
            for nlev in 72 128; do
                for nprob in 32 48 64 96 128 136 192 272 2048; do
                    for nw in 1 2 4 8 16; do
                        $standalone-$machine -np $nprob -nr $nlev -nw $nw -nc $nrhs -m cr_amxm
                        $standalone-$machine -np $nprob -nr $nlev -nw $nw -nc $nrhs -m thomas
                        $scream-$machine -np $nprob -nr $nlev -nc $nrhs -nw $nw -m cr
                        $scream-$machine -np $nprob -nr $nlev -nc $nrhs -nw $nw -m thomas
                    done
                done
            done
        done
    done
else
    for trial in $(seq $nrep); do
        echo ">>> fast"
        for nrhs in 1; do
            for nlev in 72 128; do
                for nprob in 96 136 192 272 384 768 1632 2176 4096 8192; do
                    $prefix $standalone-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas
                    $prefix $standalone-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas_pack_a1xm
                    $prefix $scream-$machine -np $nprob -nr $nlev -nc $nrhs -m thomas -1a -nop
                    $prefix $scream-$machine -np $nprob -nr $nlev -nc $nrhs -m thomas -1a
                done
            done
        done
        # Multiple RHS
        echo ">>> nrhs"
        for nrhs in 1 2 13 16 43; do
            for nlev in 72 128; do
                for nprob in 96 136 192 272 384 768 1632 2176; do
                    $prefix $standalone-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas
                    $prefix $standalone-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas_pack_a1xm
                    $prefix $scream-$machine -np $nprob -nr $nlev -nc $nrhs -m thomas -1a -nop
                    $prefix $scream-$machine -np $nprob -nr $nlev -nc $nrhs -m thomas -1a
                done
            done
        done
        # Homme problem
        echo ">>> homme"
        for nrhs in 4 8 10 16 20; do
            for nlev in 72 128; do
                for nprob in 48 96 136 272 1024 2048; do
                    $prefix $standalone-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas_amxm
                    $prefix $standalone-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas_pack_amxm
                    $prefix $scream-$machine -np $nprob -nr $nlev -nc $nrhs -m thomas -nop
                    $prefix $scream-$machine -np $nprob -nr $nlev -nc $nrhs -m thomas
                done
            done
        done
    done
fi
