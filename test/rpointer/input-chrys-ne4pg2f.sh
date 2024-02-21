outputbase=/lcrc/group/e3sm/ac.ambradl/scratch/chrys
repo=/home/ac.ambradl/e3sm
account=e3sm

model=e3sm
res=ne4pg2_oQU480
compset=F2010

machine=chrysalis
compiler=gnu
queue=debug #compute
nnode=2
ncore=32
jmake=64

unit=ndays
runit=ndays
restn=1

finstep=121
phase1=(false 2 1 10103 0)
phase2=(true 2 1 0 0)
phase3=(false 5 0 0 0)
