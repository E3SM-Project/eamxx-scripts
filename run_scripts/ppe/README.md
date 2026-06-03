# Description of files and workflow

# Multiresolution PPE Workflow

These scripts set up and launch nested multiresolution SCREAM PPEs on `pm-gpu`.

The current multiresolution PPE design uses appendable nested LHS JSON files:

| Resolution | Members | JSON file                          | Ensemble template                                     |
| ---------- | ------: | ---------------------------------- | ----------------------------------------------------- |
| ne32       |    1024 | `append_nested_ne32_1024_lhs.json` | `run.forENS.082019.ne32pg2.F2010-SCREAMv1.pm-gpu.sh`  |
| ne128      |     256 | `append_nested_ne128_256_lhs.json` | `run.forENS.082019.ne128pg2.F2010-SCREAMv1.pm-gpu.sh` |
| ne256      |     128 | `append_nested_ne256_128_lhs.json` | `run.forENS.082019.ne256pg2.F2010-SCREAMv1.pm-gpu.sh` |

The nested design satisfies:

```text
ne256 subset ⊂ ne128 subset ⊂ ne32 full ensemble
```

This allows direct comparison of identical parameter sets across resolutions.

---

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

* The ensemble templates expect all 19 PPE parameter environment variables to be provided by `create_ensembles.py`.
* The `--old-exe` argument should point to a valid prebuilt executable directory.
* The append-nested JSON files are the preferred multiresolution sampling files.
* Resolution-specific output YAMLs and remapping choices are defined inside the corresponding `run.forENS...` template scripts.

