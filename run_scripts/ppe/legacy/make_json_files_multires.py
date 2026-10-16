#!/usr/bin/env python3

import json
import math
from dataclasses import dataclass
from typing import Literal

import numpy as np

ScaleType = Literal["linear", "log"]

N_FULL = 1024
N_NE128 = 256
N_NE256 = 128
SEED = 42

ENFORCE_LAMBDA_ORDER = True  # note: this slightly breaks strict marginal LHS for lambda columns

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

def make_nested_rank_matrix(n_full, n_ne128, n_ne256, n_dim, seed):
    """
    Constructs ranks so that:
      rows 0:128 are LHS at 128 level
      rows 0:256 are LHS at 256 level
      rows 0:1024 are LHS at 1024 level

    For each dimension independently:
      - ne256 subset gets one rank in each block of 8 full strata
      - ne128 subset gets one rank in each block of 4 full strata
      - full set gets one rank in each full stratum
    """
    rng = np.random.default_rng(seed)
    ranks = np.empty((n_full, n_dim), dtype=int)

    block128 = n_full // n_ne256  # 8
    block256 = n_full // n_ne128  # 4

    for j in range(n_dim):
        used = set()

        # Assign ne256 rows: one per 8-stratum block.
        coarse128 = rng.permutation(n_ne256)
        chosen_q256 = []

        for row, c in enumerate(coarse128):
            offset = rng.integers(0, block128)
            r = c * block128 + offset
            ranks[row, j] = r
            used.add(r)

            q = r // block256
            chosen_q256.append(q)

        # Assign additional ne128 rows: complete one per 4-stratum block.
        remaining_q_by_block8 = []
        for c, q_used in zip(coarse128, chosen_q256):
            q0 = 2 * c
            q1 = 2 * c + 1
            q_other = q1 if q_used == q0 else q0
            remaining_q_by_block8.append(q_other)

        rng.shuffle(remaining_q_by_block8)

        for k, q in enumerate(remaining_q_by_block8):
            row = n_ne256 + k
            candidates = list(range(q * block256, (q + 1) * block256))
            candidates = [r for r in candidates if r not in used]
            r = rng.choice(candidates)
            ranks[row, j] = r
            used.add(r)

        # Fill remaining full rows with all unused ranks.
        remaining_ranks = np.array([r for r in range(n_full) if r not in used])
        rng.shuffle(remaining_ranks)
        ranks[n_ne128:, j] = remaining_ranks

    return ranks

def ranks_to_unit_samples(ranks, seed):
    rng = np.random.default_rng(seed + 999)
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

def validate_nested_lhs(ranks):
    checks = []

    for n_subset, block in [(N_NE256, 8), (N_NE128, 4), (N_FULL, 1)]:
        subset = ranks[:n_subset, :]
        ok_all = True

        for j in range(subset.shape[1]):
            bins = subset[:, j] // block
            ok = len(np.unique(bins)) == n_subset
            ok_all = ok_all and ok

        checks.append((n_subset, ok_all))

    return checks

def write_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f)

def main():
    n_dim = len(PARAMETERS)

    ranks = make_nested_rank_matrix(
        n_full=N_FULL,
        n_ne128=N_NE128,
        n_ne256=N_NE256,
        n_dim=n_dim,
        seed=SEED,
    )

    print("Nested LHS validation:")
    for n_subset, ok in validate_nested_lhs(ranks):
        print(f"  first {n_subset:4d} rows are LHS: {ok}")

    unit = ranks_to_unit_samples(ranks, seed=SEED)
    samples = transform_samples(unit)

    ne256_samples = samples[:N_NE256]
    ne128_samples = samples[:N_NE128]
    ne32_samples = samples[:N_FULL]

    write_json("ne32_1024_lhs.json", ne32_samples)
    write_json("ne128_256_lhs_subset.json", ne128_samples)
    write_json("ne256_128_lhs_subset.json", ne256_samples)

    metadata = {
        "seed": SEED,
        "n_full": N_FULL,
        "n_ne128": N_NE128,
        "n_ne256": N_NE256,
        "nesting": {
            "ne256_indices_in_ne32": list(range(N_NE256)),
            "ne128_indices_in_ne32": list(range(N_NE128)),
            "ne256_is_subset_of_ne128": True,
            "ne128_is_subset_of_ne32": True,
        },
        "parameter_order": [p.name for p in PARAMETERS],
        "parameter_scales": [p.scale for p in PARAMETERS],
        "enforce_lambda_low_le_lambda_high": ENFORCE_LAMBDA_ORDER,
    }

    write_json("nested_indices.json", metadata)

    print("\nWrote:")
    print("  ne32_1024_lhs.json")
    print("  ne128_256_lhs_subset.json")
    print("  ne256_128_lhs_subset.json")
    print("  nested_indices.json")

    print("\nParameter order:")
    for i, p in enumerate(PARAMETERS, start=1):
        print(f"{i:2d}. {p.name:38s} {p.scale:6s} min={p.minimum:g} max={p.maximum:g}")

if __name__ == "__main__":
    main()
