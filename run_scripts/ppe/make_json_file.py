#!/usr/bin/env python3

import json
import math
from dataclasses import dataclass
from typing import List, Literal

import numpy as np

try:
    from scipy.stats import qmc
except ImportError as e:
    raise SystemExit(
        "This script requires scipy. Install it with:\n"
        "  pip install scipy"
    ) from e


ScaleType = Literal["linear", "log"]


@dataclass(frozen=True)
class Parameter:
    number: int
    lhs_name: str
    description: str
    default: float
    minimum: float
    maximum: float
    scale: ScaleType


PARAMETERS: List[Parameter] = [
    Parameter(1,  "thl2tune",                  "temperature variance scaling factor",                                 1.0,      0.1,       10.0,      "linear"),
    Parameter(2,  "qw2tune",                   "moisture variance scaling factor",                                    1.0,      0.1,       10.0,      "linear"),
    Parameter(3,  "length_fac",                "turbulent length scale",                                              0.5,      0.1,       2.0,      "linear"),
    Parameter(4,  "c_diag_3rd_mom",            "coef for analytical formulation of 3rd moment of vertical velocity",  7.0,      0.01,      10.0,      "log"),
    Parameter(5,  "Ckh",                       "eddy diffusivity scaling factor for heat and moisture",               0.1,      0.01,      0.2,       "log"),
    Parameter(6,  "Ckm",                       "eddy diffusivity scaling factor for momentum",                        0.1,      0.01,      0.2,       "log"),
    Parameter(7,  "lambda_low",                "minimum value of stability correction",                               0.001,    0.0001,    0.1,       "linear"),
    Parameter(8,  "lambda_high",               "maximum value of stability correction",                               0.08,     0.0001,    0.1,       "linear"),
    Parameter(9,  "spa_to_nc",                 "Scaling factor for turning CCN into nc in SPA",                       2000.0,   100.0,     4000.0,    "linear"),
    Parameter(10, "eci",                       "liquid/ice collision/collection coefficient",                         0.5,      0.1,       1.0,       "linear"),
    Parameter(11, "eri",                       "ice/rain collision/collection coefficient",                           1.0,      0.1,       1.0,       "linear"),
    Parameter(12, "k_acc",                     "scaling factor on accretion",                                         67.0,     0.01,      100.0,     "log"),
    Parameter(13, "dep_nuc_exponent",          "deposition nucleation exponent",                                      0.304,    0.2,       0.304,     "linear"),
    Parameter(14, "max_total_ni",              "limiter on max ni value in a cell",                                   7.4e5,    5.0e5,     2.0e6,     "log"),
    Parameter(15, "ice_sed_knob",              "ice fall speed",                                                      1.0,      0.8,       2.0,       "linear"),
    Parameter(16, "D_breakup_cutoff",          "rain self collection and breakup",                                    280e-6,   0.0,       500e-6,    "linear"),
    Parameter(17, "autoconversion_prefactor",  "autoconversion prefactor",                                            2700.0,   10.0,      20000.0,   "log"),
    Parameter(18, "autoconversion_qc_exponent","autoconversion qc exponent",                                          2.47,     2.0,       4.0,       "linear"),
    Parameter(19, "autoconversion_radius",     "autoconversion radius",                                               25e-6,    25e-6,     50e-6,     "linear"),
]


def scale_unit_to_range(u: float, p: Parameter) -> float:
    """
    Map a value u in [0,1) to the parameter's physical range.
    """
    if not (0.0 <= u <= 1.0):
        raise ValueError(f"u must be in [0,1], got {u}")

    if p.scale == "linear":
        return p.minimum + u * (p.maximum - p.minimum)

    if p.scale == "log":
        if p.minimum <= 0 or p.maximum <= 0:
            raise ValueError(f"Log-scaled parameter must have positive bounds: {p.lhs_name}")
        log_min = math.log10(p.minimum)
        log_max = math.log10(p.maximum)
        return 10 ** (log_min + u * (log_max - log_min))

    raise ValueError(f"Unknown scale type {p.scale!r} for parameter {p.lhs_name}")


def enforce_lambda_order(row: List[float]) -> List[float]:
    """
    Ensure lambda_low <= lambda_high within each sampled parameter set.

    Indices are 0-based:
      lambda_low  -> column 7  -> index 6
      lambda_high -> column 8  -> index 7
    """
    lambda_low_idx = 6
    lambda_high_idx = 7

    if row[lambda_low_idx] > row[lambda_high_idx]:
        row[lambda_low_idx], row[lambda_high_idx] = row[lambda_high_idx], row[lambda_low_idx]

    return row


def generate_lhs_samples(
    n_samples: int,
    seed: int = 42,
    scramble: bool = True,
) -> List[List[float]]:
    """
    Generate Latin hypercube samples in the exact parameter order above.
    """
    n_dim = len(PARAMETERS)
    sampler = qmc.LatinHypercube(d=n_dim, scramble=scramble, seed=seed)
    unit_samples = sampler.random(n=n_samples)

    scaled_samples: List[List[float]] = []
    for row in unit_samples:
        scaled_row = [
            scale_unit_to_range(float(u), p)
            for u, p in zip(row, PARAMETERS)
        ]

        scaled_row = enforce_lambda_order(scaled_row)
        scaled_samples.append(scaled_row)

    return scaled_samples


def main() -> None:
    n_samples = 256
    seed = 42
    out_file = "normranked_LH_sampling_base10_updated_reduced_test.json"

    samples = generate_lhs_samples(n_samples=n_samples, seed=seed, scramble=True)

    with open(out_file, "w") as f:
        json.dump(samples, f)

    print(f"Wrote {len(samples)} samples to {out_file}")
    print("Parameter order:")
    for p in PARAMETERS:
        print(f"{p.number:2d}. {p.lhs_name:28s} [{p.scale}]  min={p.minimum} max={p.maximum}")

    # sanity check
    n_bad = sum(1 for row in samples if row[6] > row[7])
    print(f"\nNumber of sets with lambda_low > lambda_high after correction: {n_bad}")


if __name__ == "__main__":
    main()
