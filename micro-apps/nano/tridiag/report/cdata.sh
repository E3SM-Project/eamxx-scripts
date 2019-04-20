machine=$1 # v100, skx, knl
nrep=101 # extreme; went with >= 25

prefix=""
if [ $machine == "knl" ]; then
    export OMP_NUM_THREADS=136 # or 272
    export OMP_PROC_BIND=false
    export KMP_AFFINITY=balanced,granularity=fine
    prefix=" numactl -i 1"
elif [ $machine == "v100" ]; then
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
        for nlev in 72 128 256; do
            for nprob in 96 136 192 272 384 768 1632 2176 2400 2800 3200 4800 6400 12800 25600; do
                ./tridiag -np $nprob -nr $nlev -nc 1 -m cusparse
                for nw in 1 2 4 8; do
                    ./tridiag -np $nprob -nr $nlev -nw $nw -nc 1 -m cr_a1x1
                    ./tridiag -np $nprob -nr $nlev -nw $nw -nc 1 -m cr_a1x1p
                    ./tridiag -np $nprob -nr $nlev -nw $nw -nc 1 -m cr_a1xm
                    ./tridiag -np $nprob -nr $nlev -nw $nw -nc 1 -m cr_amxm
                    ./tridiag -np $nprob -nr $nlev -nw $nw -nc 1 -m thomas
                done
            done
        done
        # Multiple RHS
        echo ">>> nrhs"
        for nrhs in 2 3 4 8 13 16 32 43 64; do
            for nlev in 72 128 256; do
                for nprob in 96 136 192 272 384 768 1632 2176 2400 2800 3200 4800 6400 12800 25600; do
                    ./tridiag -np $nprob -nr $nlev -nc $nrhs -m cusparse
                    for nw in 1 2 4 8 16; do
                        ./tridiag -np $nprob -nr $nlev -nw $nw -nc $nrhs -m cr_a1xm
                        ./tridiag -np $nprob -nr $nlev -nw $nw -nc $nrhs -m thomas
                    done
                done
            done
        done
        # Homme problem
        echo ">>> homme"
        for nrhs in 4 8 10 16 20; do
            for nlev in 72 128 256; do
                for nprob in 32 48 64 96 128 136 192 256 272 512 800 1024 1600 2048; do
                    ./tridiag -np $nprob -nr $nlev -nc $nrhs -m cusparse
                    for nw in 1 2 4 8 16; do
                        ./tridiag -np $nprob -nr $nlev -nw $nw -nc $nrhs -m cr_amxm
                        ./tridiag -np $nprob -nr $nlev -nw $nw -nc $nrhs -m thomas
                    done
                done
            done
        done
    done
else
    for trial in $(seq $nrep); do
        # Multiple RHS
        echo ">>> nrhs"
        for nrhs in 1 2 3 4 8 13 16 32 43 64; do
            for nlev in 72 128 256; do
                for nprob in 96 136 192 272 384 768 1632 2176 2400 2800 3200 4800 6400 12800 25600; do
                    $prefix ./tridiag-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m gttr
                    $prefix ./tridiag-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m dttr
                    $prefix ./tridiag-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas
                    $prefix ./tridiag-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas_pack_a1xm
                done
            done
        done
        # Homme problem
        echo ">>> homme"
        for nrhs in 4 8 10 16 20; do
            for nlev in 72 128 256; do
                for nprob in 32 48 64 96 128 136 192 256 272 512 800 1024 1600 2048; do
                    $prefix ./tridiag-$machine -ng -np $(expr $nrhs \* $nprob) -nr $nlev -nc 1 -m gttr
                    $prefix ./tridiag-$machine -ng -np $(expr $nrhs \* $nprob) -nr $nlev -nc 1 -m dttr
                    $prefix ./tridiag-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas_amxm
                    $prefix ./tridiag-$machine -ng -np $nprob -nr $nlev -nc $nrhs -m thomas_pack_amxm
                done
            done
        done
    done
fi
