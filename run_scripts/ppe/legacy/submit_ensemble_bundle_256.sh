#!/bin/csh                                                                                                                                                                                                                                                                     

#SBATCH --account=e3sm                       
#SBATCH -C gpu
#SBATCH --job-name=bundletryall
#SBATCH -q regular                                                  
#SBATCH -N 256
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=none
#SBATCH -t 19:29:00    
#SBATCH -o sout.%j     
#SBATCH -e serr.%j
#SBATCH --mail-type=BEGIN                      
#SBATCH --mail-type=END                        
###SBATCH --mail-user=USEREMAIL

echo "SLURM_JOB_ID="$SLURM_JOB_ID
echo "SLURM_JOB_NAME="$SLURM_JOB_NAME
echo "SLURM_LOCALID="$SLURM_LOCALID
echo "SLURM_NODEID="$SLURM_NODEID

pwd

set base = "/pscratch/sd/t/terai/e3sm_scratch/pm-gpu"
set caseprefix = "PPEensemble_1node_full256.ne32pg2_ne32pg2.F2010-SCREAMv1.20260319"

foreach i (`seq 0 255`)
  set mem = `printf "m%03d" $i`
  set dir = "${base}/${caseprefix}.${mem}/case_scripts"

  date
  cd $dir
  ./case.submit --no-batch -v >& bsubmitout.txt &
end

wait
echo "Done!"
date
