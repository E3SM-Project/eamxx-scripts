outputbase=/lcrc/group/e3sm/ac.ambradl/scratch/chrys
repo=/home/ac.ambradl/SCREAMref
account=e3sm

model=scream
res=ne30pg2_ne30pg2
compset=F2010-SCREAMv1

machine=chrysalis
compiler=gnu
queue=compute
nnode=5
ncore=64
jmake=64

unit=nhours
runit=nhours
restn=6

finstep=95
phase1=(false 24 0 010101 64800)
phase2=(true 12 2 0 0)
phase3=(false 48 0 0 0)
