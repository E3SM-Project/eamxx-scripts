outputbase=/lustre/orion/cli115/proj-shared/ambradl/e3sm_scratch/
repo=/ccs/home/ambradl/repo/SCREAMref
account=cli115

model=scream
res=ne4pg2_ne4pg2
compset=F2010-SCREAMv1

machine=frontier-scream-gpu
compiler=crayclang-scream
queue=batch
nnode=1
ncore=8
jmake=8

unit=ndays
runit=ndays
restn=1

finstep=119
phase1=(false 2 1 10103 0)
phase2=(true 2 1 0 0)
phase3=(false 5 0 0 0)
