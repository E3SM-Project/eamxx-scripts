outputbase=/lustre/orion/cli115/proj-shared/ambradl/e3sm_scratch/
repo=/ccs/home/ambradl/repo/SCREAMref
account=cli115

model=scream
res=ne30pg2_ne30pg2
compset=F2010-SCREAMv1

machine=frontier-scream-gpu
compiler=crayclang-scream
queue=batch
nnode=1
ncore=8
jmake=8

unit=nhours
runit=nhours
restn=6

finstep=95
phase1=(false 24 0 010101 64800)
phase2=(true 12 2 0 0)
phase3=(false 48 0 0 0)
