# Multiresolution PPE Workflow

These scripts set up and launch nested multiresolution SCREAM PPEs on `pm-gpu`.

The PPE design uses appendable nested LHS JSON files:

| Resolution | Members | Nodes per member | JSON file                          |
| ---------- | ------: | ---------------: | ---------------------------------- |
| ne32       |    1024 |                1 | `append_nested_ne32_1024_lhs.json` |
| ne128      |     256 |               16 | `append_nested_ne128_256_lhs.json` |
| ne256      |     128 |               64 | `append_nested_ne256_128_lhs.json` |

The nested design satisfies:

```text
ne256 subset ⊂ ne128 subset ⊂ ne32 full ensemble
```

This allows direct comparison of identical parameter sets across resolutions.

There are three sets of scripts in this directory:

* **Unified workflow (current):** one run script, `run.unified.F2010-SCREAMv1.pm-gpu.sh`, for every resolution and for both the control and the ensemble members, used with `create_ensembles_unified.py`, `make_json_files_multires_unified.py` and the shared parameter list `ppe_parameters.py`. Any JSON file works: SHOC + P3, ZM only, or combined.
* **First ZM scripts:** `run.ZM.082019.F2010-SCREAMv1.pm-gpu.sh` with `create_ensembles.py` and `make_json_files_multires_append.py`. They vary the fixed 19 SHOC + P3 parameters only, and are kept until the unified scripts replace them.
* **Previous workflow:** separate `run.default...` and `run.forENS...` scripts per resolution, kept as the record of the earlier PPE. See [Previous workflow](#previous-workflow-per-resolution-scripts) at the end of this file.

---

# Unified workflow

`run.unified.F2010-SCREAMv1.pm-gpu.sh` (ZM + SPC, L128v4, Aug 2019 start) covers ne32, ne128 and ne256, as either the control or an ensemble member:

* **Resolution** comes from the first argument (`ne32`, `ne128` or `ne256`) or from `RES` in the environment. `create_ensembles_unified.py` passes no arguments, so set `RES` when calling it.
* **Mode:** if `MEMBER_ID` is set, which `create_ensembles_unified.py` always does, the script creates an ensemble member. Otherwise it creates the control.
* Everything that differs between resolutions (PE layout, walltime, SPC remap, output remap, land initial conditions, `ATM_NCPL`) is in `resolution_settings()` in the script.

Case locations, set by `CONTROL_ROOT_BASE` in the script:

```text
control: $SCRATCH/e3sm_scratch/pm-gpu/<res>_ppe_v4_config_new_workflow_unified/PPEensemble_<N>node.<res>pg2_<res>pg2.F2010-SCREAMv1.<CHECKOUT>
members: $SCRATCH/e3sm_scratch/pm-gpu/<res>_ppe_v4_config_new_workflow_unified_ens/<control case name>.mXXX
```

## 1. Choose the parameters and make the JSON files

Every parameter the PPE can vary is listed in `ppe_parameters.py`, with its namelist key, range, scale, kind and code default. `python ppe_parameters.py` prints the list. Which parameters an ensemble actually varies is decided by its JSON file: anything not in the file runs at its code default.

Generate a nested set of JSON files (ne256 ⊂ ne128 ⊂ ne32) for any selection of schemes and parameters:

```bash
python make_json_files_multires_unified.py --schemes zm --tag zm_v1
python make_json_files_multires_unified.py --schemes shoc,p3,zm --exclude p3::max_total_ni,shoc::c_diag_3rd_mom --tag combined_v1
```

This writes `append_nested_ne32_1024_lhs_<tag>.json`, `append_nested_ne128_256_lhs_<tag>.json`, `append_nested_ne256_128_lhs_<tag>.json` and `append_nested_lhs_metadata_<tag>.json`. Existing files are not overwritten unless you pass `--force`.

The files carry their parameter names:

```json
{"parameters": ["zm::tau", "zm::alfa", "..."], "samples": [[...], [...]]}
```

Integer parameters, such as `zm::num_cin`, are sampled as equal-probability integer levels between their minimum and maximum.

The original 19-column files (SHOC + P3, plain lists of values) still work everywhere. `--schemes shoc,p3 --tag v2 --legacy-format --force` reproduces the `_v2` files exactly.

To make a new parameter available, add one line to `ppe_parameters.py`. No other script needs to change.

## 2. Run the control

The control builds the executable that every member of that resolution reuses.

```bash
./run.unified.F2010-SCREAMv1.pm-gpu.sh ne128
```

By default the control builds and submits to `regular` with the walltime from `resolution_settings()`. Useful overrides:

```bash
Q=debug WALLTIME=00:30:00 ./run.unified.F2010-SCREAMv1.pm-gpu.sh ne32   # quick test
do_case_build=false ./run.unified.F2010-SCREAMv1.pm-gpu.sh ne32         # reuse this control's existing build
```

Don't re-run a control while its job is active: it edits and resubmits the same case.

## 3. Create member cases

Create members with `create_ensembles_unified.py`, with `--old-exe` pointing at the control's `build` directory:

```bash
RES=ne128 python create_ensembles_unified.py \
  --template ./run.unified.F2010-SCREAMv1.pm-gpu.sh \
  --json ./append_nested_ne128_256_lhs_zm_v1.json \
  --old-exe <ne128 control case>/build \
  --start 0 \
  --count 32
```

`--start` is the first member index and `--count` is the number of members. Member IDs are zero-indexed, so `--start 0 --count 32` gives `m000` through `m031`.

Members are created but not submitted. To submit each member as its own job instead of bundling, put `do_case_submit=true` in front of the command.

To hold some parameters in the file at their defaults, or to check a file first:

```bash
--exclude-scheme p3              # hold every p3:: parameter at its default
--exclude zm::mcsp_q_coeff       # hold one parameter at its default (repeatable, or comma-separated)
--dry-run                        # check the file and print each member's changes; creates nothing
```

Before creating anything, `create_ensembles_unified.py` checks the selected rows:

* every parameter name is in `ppe_parameters.py`;
* every row has the right length;
* integer parameters are whole numbers;
* values are inside the ranges in `ppe_parameters.py`. This is an error for named files unless you pass `--allow-out-of-range`, and a warning for the original 19-column files, since earlier designs used wider ranges.

It then prints which parameters vary and which are held at their defaults. Each member receives its changes as one string, `PPE_CHANGES="key=value ..."`.

Member creation stops with an error if `RES` is missing, if `--old-exe` has no `e3sm.exe`, if `--old-exe` is not a build for that resolution, or if `atmchange` rejects a PPE parameter name. Each member's log (`ensemble_logs/mXXX.log`) ends with the process list and the parameter values actually in its namelist, and each case records its changes in `case_scripts/PPE_PARAMETERS.txt`.

### Creating members in parallel

Each member takes about a minute to create. Members are independent, so run several `create_ensembles_unified.py` processes on non-overlapping ranges. For example, all 1024 ne32 members in 32 chunks:

```bash
export RES=ne32
B=<ne32 control case>/build
N=32                                # chunks running at once
STEP=$(( (1024 + N - 1) / N ))      # members per chunk
mkdir -p chunk_logs
for s in $(seq 0 $STEP 1023); do
  python create_ensembles_unified.py --template ./run.unified.F2010-SCREAMv1.pm-gpu.sh \
    --json ./append_nested_ne32_1024_lhs.json --old-exe $B \
    --start $s --count $STEP > chunk_logs/chunk_$s.log 2>&1 &
done
wait
```

On a login node, keep `N` around 8. For more, run the loop on an interactive CPU node, for example `salloc -N 1 -C cpu -q interactive -t 01:00:00 -A e3sm`. Make sure `python` is available in that shell (`module load python`).

### Checking progress

```bash
grep -l "All done" $(grep -l "MODE=ensemble" ensemble_logs/m*.log) | wc -l   # members finished
tail -qn1 chunk_logs/chunk_*.log                                              # member each chunk is on
grep -l RuntimeError chunk_logs/*.log                                         # chunks that stopped on a failure
grep -L "All done" ensemble_logs/m*.log                                       # members that did not finish
```

### Notes on member creation

* `create_ensembles_unified.py` writes `ensemble_logs/<member>.log` into the directory it is run from, named by member ID only. Use a separate working directory per resolution, or move `ensemble_logs/` before starting the next resolution, so logs are not overwritten.
* To redo a member, delete its case directory first. CIME's `create_newcase` refuses an existing case directory.
* A chunk stops at its first failed member while the other chunks keep going. Fix the problem, delete the failed member's case directory, and rerun that range.

## 4. Submit members

Submit members in bundles with a csh batch script that runs `./case.submit --no-batch` in each member's `case_scripts` directory inside one allocation. `submit_ensemble_ne32_ZM.sh` submits ne32 `m000` through `m127` on 128 nodes:

```bash
sbatch submit_ensemble_ne32_ZM.sh
```

For other resolutions or member ranges, update:

* `#SBATCH -N`: nodes per member × members in the bundle
* walltime
* `base`: the members directory (`<res>_ppe_v4_config_new_workflow_unified_ens`)
* `caseprefix`: the control case name
* the member range in the `foreach` loop

Before a large bundle, run one member on debug:

```bash
cd <members directory>/<control case name>.m000/case_scripts
./case.submit -a="-t 00:30:00 --qos=debug"
```

## Notes

* The run script has no parameter list of its own. It applies whatever `create_ensembles_unified.py` passes in `PPE_CHANGES`.
* The run script changes `physics::atm_procs_list` without `-b`, so the namelist contains `zm` and `spc` right away. Otherwise `zm::` parameters can't be checked when they're applied.
* At this code checkout `spc_remap_file` has no default for ne32 or ne128, so the script sets it explicitly for every resolution.
* All resolutions use `set_cld_frac_r_to_one=True`. The previous ne256 scripts did not set it.

---

# Previous workflow (per-resolution scripts)

The earlier PPE used one control script (`run.default...`) and one ensemble template (`run.forENS...`) per resolution:

| Resolution | Ensemble template                                     |
| ---------- | ----------------------------------------------------- |
| ne32       | `run.forENS.082019.ne32pg2.F2010-SCREAMv1.pm-gpu.sh`  |
| ne128      | `run.forENS.082019.ne128pg2.F2010-SCREAMv1.pm-gpu.sh` |
| ne256      | `run.forENS.082019.ne256pg2.F2010-SCREAMv1.pm-gpu.sh` |

## 1. Build or prepare a baseline executable

Each PPE member reuses an existing executable through the `--old-exe` argument to `create_ensembles.py`.

Build or identify a baseline case for the desired resolution, then use its `build` directory as the `--old-exe` path.

Example ne32 executable:

```bash
/pscratch/sd/b/beydoun/e3sm_scratch/pm-gpu/ne32_ppe_prod/PPEensemble_1node.ne32pg2_ne32pg2.F2010-SCREAMv1.20260515/build
```

Example ne128 executable:

```bash
/pscratch/sd/b/beydoun/e3sm_scratch/pm-gpu/ne128_ppe_prod/PPEensemble_16node.ne128pg2_ne128pg2.F2010-SCREAMv1.20260515/build
```

Example ne256 executable:

```bash
/pscratch/sd/b/beydoun/e3sm_scratch/pm-gpu/ne256_ppe_prod/PPEensemble_16node.ne256pg2_ne256pg2.F2010-SCREAMv1.20260515/build
```

Adjust these paths as needed.

---

## 2. Generate PPE member cases

The script `create_ensembles.py` reads:

* a JSON file containing the PPE parameter values
* a resolution-specific ensemble template
* an existing executable path
* an optional `--start` and `--count` range

`--start` is the first member index.
`--count` is the number of members to create.

The generated member IDs are zero-indexed:

```text
--start 0 --count 512    -> m000 through m511
--start 512 --count 512  -> m512 through m1023
```

---

## 3. ne32: create 1024 member cases

### First half

```bash
python create_ensembles.py \
  --json ./append_nested_ne32_1024_lhs.json \
  --template ./run.forENS.082019.ne32pg2.F2010-SCREAMv1.pm-gpu.sh \
  --old-exe /pscratch/sd/b/beydoun/e3sm_scratch/pm-gpu/ne32_ppe_prod/PPEensemble_1node.ne32pg2_ne32pg2.F2010-SCREAMv1.20260515/build \
  --start 0 \
  --count 512
```

### Second half

```bash
python create_ensembles.py \
  --json ./append_nested_ne32_1024_lhs.json \
  --template ./run.forENS.082019.ne32pg2.F2010-SCREAMv1.pm-gpu.sh \
  --old-exe /pscratch/sd/b/beydoun/e3sm_scratch/pm-gpu/ne32_ppe_prod/PPEensemble_1node.ne32pg2_ne32pg2.F2010-SCREAMv1.20260515/build \
  --start 512 \
  --count 512
```

---

## 4. ne128: create 256 member cases

Example 32-member batch:

```bash
python create_ensembles.py \
  --json ./append_nested_ne128_256_lhs.json \
  --template ./run.forENS.082019.ne128pg2.F2010-SCREAMv1.pm-gpu.sh \
  --old-exe /pscratch/sd/b/beydoun/e3sm_scratch/pm-gpu/ne128_ppe_prod/PPEensemble_16node.ne128pg2_ne128pg2.F2010-SCREAMv1.20260515/build \
  --start 0 \
  --count 32
```

Repeat with:

```text
--start 0    --count 32
--start 32   --count 32
--start 64   --count 32
--start 96   --count 32
--start 128  --count 32
--start 160  --count 32
--start 192  --count 32
--start 224  --count 32
```

---

## 5. ne256: create 128 member cases

Example 16-member batch:

```bash
python create_ensembles.py \
  --json ./append_nested_ne256_128_lhs.json \
  --template ./run.forENS.082019.ne256pg2.F2010-SCREAMv1.pm-gpu.sh \
  --old-exe /pscratch/sd/b/beydoun/e3sm_scratch/pm-gpu/ne256_ppe_prod/PPEensemble_16node.ne256pg2_ne256pg2.F2010-SCREAMv1.20260515/build \
  --start 0 \
  --count 16
```

Repeat as needed with:

```text
--start 0    --count 16
--start 16   --count 16
--start 32   --count 16
--start 48   --count 16
--start 64   --count 16
--start 80   --count 16
--start 96   --count 16
--start 112  --count 16
```

If using a different bundle size, adjust `--count` accordingly.

---

## 6. Submit ensemble bundles

After cases are created, use a bundle submission script such as:

```bash
sbatch submit_ensemble_bundle_256.sh
```

Before submitting, update the bundle script for the target resolution:

* `#SBATCH -N`
* walltime
* `base`
* `caseprefix`
* member index range
* email address, if desired

The number of allocated nodes should match the number of concurrent member cases submitted when each member uses one node. For multi-node cases, adjust the node count and member range accordingly.

---

## 7. Notes

* The ensemble templates expect all 19 PPE parameter environment variables to be provided by `create_ensembles.py`. `create_ensembles_unified.py` also sets them for the original 19-column files, so these templates work with it too.
* The `--old-exe` argument should point to a valid prebuilt executable directory.
* The append-nested JSON files are the preferred multiresolution sampling files.
* Resolution-specific output YAMLs and remapping choices are defined inside the corresponding `run.forENS...` template scripts.
