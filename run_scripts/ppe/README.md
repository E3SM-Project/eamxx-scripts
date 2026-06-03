# Description of files and workflow

Note that these were performed on pm-gpu and might need modifications for running elsewhere. The most up to date json files are the append_nested_*_lhs.json files generated using python make_json_files_multires_append.py

1) Build (and test) baseline case to create the e3sm.exe that all PPE ensemble members will use (`run.default.ne32pg2.F2010-SCREAMv1.pm-gpu.sh`). This uses 082019 initial conditions and designed to run for 13 months. 

2) Set the simulation length, outputs, etc in the template script that will create each of the ensemble cases (`run.forENS.082019.ne32pg2.F2010-SCREAMv1.pm-gpu.sh`).

3) Create all PPE member cases (`create_ensembles.py`)
When running the python script, provide it with flags pointing to the json file with the tuning parameters, the template script to make each of the ensemble cases, and the location of the e3sm.exe.

```bash
python create_ensembles.py --json ./append_nested_ne32_1024_lhs.json --template ./run.forENS.082019.ne32pg2.F2010-SCREAMv1.pm-gpu.sh  --old-exe /pscratch/sd/b/beydoun/e3sm_scratch/pm-gpu/ne32_ppe_prod/PPEensemble_1node.ne32pg2_ne32pg2.F2010-SCREAMv1.20260515/build --start 14 --count 12
```
Note that --start and --count indicate what ppe member to start with and how many to process. 

4) Modify the case directory, case name, GPU counts, job time, username, and in the bundle script (`submit_ensemble_bundle_256.sh`) and submit.

```bash
sbatch submit_ensemble_bundle_256.sh
```
