#!/usr/bin/env python3

import json
import math
from dataclasses import dataclass
from typing import Literal

import numpy as np

ScaleType = Literal["linear", "log"]

N_NE256 = 128
N_NE128 = 256
N_NE32 = 1024
SEED = 42

ENFORCE_LAMBDA_ORDER = True

OUT_NE256 = "append_nested_ne256_128_lhs.json"
OUT_NE128 = "append_nested_ne128_256_lhs.json"
OUT_NE32 = "append_nested_ne32_1024_lhs.json"
OUT_META = "append_nested_lhs_metadata.json"


@dataclass(frozen=True)
class Parameter:
    name: str
    minimum: float
    maximum: float
    scale: ScaleType


PARAMETERS = [
    Parameter("thl2tune", 0.1, 10.0, "linear"),
    Parameter("qw2tune", 0.1, 10.0, "linear"),
    Parameter("length_fac", 0.1, 2.0, "linear"),
    Parameter("c_diag_3rd_mom", 0.01, 10.0, "log"),
    Parameter("coeff_kh", 0.01, 0.2, "log"),
    Parameter("coeff_km", 0.01, 0.2, "log"),
    Parameter("lambda_low", 0.0001, 0.1, "linear"),
    Parameter("lambda_high", 0.0001, 0.1, "linear"),
    Parameter("spa_ccn_to_nc_factor", 100.0, 4000.0, "linear"),
    Parameter("cldliq_to_ice_collection_factor", 0.1, 1.0, "linear"),
    Parameter("rain_to_ice_collection_factor", 0.1, 1.0, "linear"),
    Parameter("accretion_prefactor", 0.01, 100.0, "log"),
    Parameter("deposition_nucleation_exponent", 0.2, 0.304, "linear"),
    Parameter("max_total_ni", 5.0e5, 1.8e6, "log"),
    Parameter("ice_sedimentation_factor", 0.7, 2.0, "linear"),
    Parameter("rain_selfcollection_breakup_diameter", 0.0, 500e-6, "linear"),
    Parameter("autoconversion_prefactor", 10.0, 20000.0, "log"),
    Parameter("autoconversion_qc_exponent", 2.0, 4.0, "linear"),
    Parameter("autoconversion_radius", 25e-6, 50e-6, "linear"),
]


def scale_value(u, p):
    if p.scale == "linear":
        return p.minimum + u * (p.maximum - p.minimum)

    log_min = math.log10(p.minimum)
    log_max = math.log10(p.maximum)
    return 10 ** (log_min + u * (log_max - log_min))


def build_appendable_ranks(n_full, n_mid, n_core, n_dim, seed):
    """
    Build ranks in append order:

      rows 0:128     are valid 128-point LHS
      rows 0:256     are valid 256-point LHS
      rows 0:1024    are valid 1024-point LHS

    Ranks are defined on the finest 1024-stratum grid.

    For each dimension:
      - First 128 rows: one per 8-bin block
      - Next 128 rows: fill the missing 4-bin half inside each 8-bin block
      - Next 768 rows: fill remaining 1024 fine bins
    """
    rng = np.random.default_rng(seed)
    ranks = np.empty((n_full, n_dim), dtype=int)

    fine_per_core = n_full // n_core  # 8
    fine_per_mid = n_full // n_mid    # 4

    for j in range(n_dim):
        used = set()

        # -----------------------------
        # Stage 1: core 128-point LHS
        # -----------------------------
        core_blocks = rng.permutation(n_core)

        used_mid_blocks = []

        for row, core_block in enumerate(core_blocks):
            # core block spans 8 fine strata
            start = core_block * fine_per_core
            end = start + fine_per_core

            # randomly choose one of the 8 fine strata
            r = int(rng.integers(start, end))

            ranks[row, j] = r
            used.add(r)

            # record which 256-level 4-bin block this occupied
            used_mid_blocks.append(r // fine_per_mid)

        # -----------------------------
        # Stage 2: append 128 to make 256 LHS
        # -----------------------------
        # Each 8-bin core block contains two 4-bin mid blocks.
        # The appended point fills the other 4-bin block.
        append_mid_blocks = []

        for core_block, used_mid in zip(core_blocks, used_mid_blocks):
            mid0 = 2 * core_block
            mid1 = 2 * core_block + 1
            other_mid = mid1 if used_mid == mid0 else mid0
            append_mid_blocks.append(other_mid)

        rng.shuffle(append_mid_blocks)

        for k, mid_block in enumerate(append_mid_blocks):
            row = n_core + k
            start = mid_block * fine_per_mid
            end = start + fine_per_mid

            candidates = [r for r in range(start, end) if r not in used]
            r = int(rng.choice(candidates))

            ranks[row, j] = r
            used.add(r)

        # -----------------------------
        # Stage 3: append 768 to make 1024 LHS
        # -----------------------------
        remaining = np.array([r for r in range(n_full) if r not in used])
        rng.shuffle(remaining)

        ranks[n_mid:, j] = remaining

    return ranks


def ranks_to_unit(ranks, seed):
    rng = np.random.default_rng(seed + 12345)
    jitter = rng.random(ranks.shape)
    return (ranks + jitter) / ranks.shape[0]


def transform_samples(unit):
    samples = []

    for row in unit:
        vals = [scale_value(float(u), p) for u, p in zip(row, PARAMETERS)]

        if ENFORCE_LAMBDA_ORDER:
            lo_idx = 6
            hi_idx = 7
            if vals[lo_idx] > vals[hi_idx]:
                vals[lo_idx], vals[hi_idx] = vals[hi_idx], vals[lo_idx]

        samples.append(vals)

    return samples


def validate_ranks(ranks):
    """
    Check strict nested LHS property in rank space before lambda swapping.
    """
    checks = []

    for n_subset in [N_NE256, N_NE128, N_NE32]:
        subset = ranks[:n_subset]
        block = N_NE32 // n_subset

        ok_all = True
        for j in range(subset.shape[1]):
            coarse_bins = subset[:, j] // block
            ok = len(np.unique(coarse_bins)) == n_subset
            ok_all = ok_all and ok

        checks.append((n_subset, ok_all))

    return checks


def write_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f)


def main():
    n_dim = len(PARAMETERS)

    ranks = build_appendable_ranks(
        n_full=N_NE32,
        n_mid=N_NE128,
        n_core=N_NE256,
        n_dim=n_dim,
        seed=SEED,
    )

    print("Rank-space nested LHS validation:")
    for n_subset, ok in validate_ranks(ranks):
        print(f"  first {n_subset:4d} rows are LHS: {ok}")

    unit = ranks_to_unit(ranks, SEED)
    samples = transform_samples(unit)

    ne256 = samples[:N_NE256]
    ne128 = samples[:N_NE128]
    ne32 = samples[:N_NE32]

    write_json(OUT_NE256, ne256)
    write_json(OUT_NE128, ne128)
    write_json(OUT_NE32, ne32)

    metadata = {
        "method": "appendable_nested_lhs",
        "seed": SEED,
        "n_ne256": N_NE256,
        "n_ne128": N_NE128,
        "n_ne32": N_NE32,
        "files": {
            "ne256": OUT_NE256,
            "ne128": OUT_NE128,
            "ne32": OUT_NE32,
        },
        "nesting": {
            "ne256_rows_in_ne32": [0, N_NE256 - 1],
            "ne128_rows_in_ne32": [0, N_NE128 - 1],
            "ne256_is_subset_of_ne128": True,
            "ne128_is_subset_of_ne32": True,
        },
        "parameter_order": [p.name for p in PARAMETERS],
        "parameter_scales": [p.scale for p in PARAMETERS],
        "enforce_lambda_low_le_lambda_high": ENFORCE_LAMBDA_ORDER,
        "note": (
            "Strict LHS is validated in rank space before optional lambda_low/lambda_high swapping. "
            "If lambda ordering is enforced, lambda_low and lambda_high marginals may no longer be strict LHS."
        ),
    }

    write_json(OUT_META, metadata)

    print("\nWrote:")
    print(f"  {OUT_NE256}")
    print(f"  {OUT_NE128}")
    print(f"  {OUT_NE32}")
    print(f"  {OUT_META}")

    print("\nAppend structure:")
    print("  rows 0-127   -> ne256")
    print("  rows 0-255   -> ne128")
    print("  rows 0-1023  -> ne32")

    print("\nParameter order:")
    for i, p in enumerate(PARAMETERS, start=1):
        print(f"{i:2d}. {p.name:38s} {p.scale:6s} min={p.minimum:g} max={p.maximum:g}")


if __name__ == "__main__":
    main()
