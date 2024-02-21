outputbase=/lcrc/group/e3sm/ac.ambradl/scratch/chrys
repo=/home/ac.ambradl/e3sm
account=e3sm

model=e3sm
res=ne30pg2_r05_EC30to60E2r2
compset=WCYCL1850

machine=chrysalis
compiler=gnu
queue=debug
nnode=5
ncore=64
jmake=64

unit=nhours
runit=nhours
restn=6

finstep=97
phase1=(false 24 0 010101 64800)
phase2=(true 12 2 0 0)
phase3=(false 48 0 0 0)
