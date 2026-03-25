# Description of files and workflow

Note that these were performed on pm-gpu and might need modifications for running elsewhere.

1) Build (and test) baseline case to create the e3sm.exe that all PPE ensemble members will use (`run.test.ne32pg2.F2010-SCREAMv1.pm-gpu.sh`)

2) Set the simulation length, outputs, etc in the template script that will create each of the ensemble cases (`run.forENS.ne32pg2.F2010-SCREAMv1.pm-gpu.sh`).

3) Create all PPE member cases (`create_ensembles.py`)
When running the python script, provide it with flags pointing to the json file with the tuning parameters, the template script to make each of the ensemble cases, and the location of the e3sm.exe.

```bash
python create_ensembles.py --json ./normranked_LH_sampling_base10.json --template ./run.forENS.ne32pg2.F2010-SCREAMv1.pm-gpu.sh --old-exe /pscratch/sd/t/terai/e3sm_scratch/pm-gpu/PPEensemble_1node.ne32pg2_ne32pg2.F2010-SCREAMv1.20260319/build
```

4) Modify the case directory, case name, GPU counts, job time, username, and in the bundle script (`submit_ensemble_bundle_256.sh`) and submit.

```bash
sbatch submit_ensemble_bundle_256.sh
```
