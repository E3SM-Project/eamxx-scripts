#!/usr/bin/env python3
"""
Appendable nested LHS for the multiresolution PPE.

Samples any selection of the parameters in ppe_parameters.py and writes one JSON file per
resolution plus a metadata file:

  rows 0-127   -> ne256      rows 0-255 -> ne128      rows 0-1023 -> ne32

Examples:
  python make_json_files_multires_append.py --schemes zm --tag zm_v1
  python make_json_files_multires_append.py --schemes shoc,p3,zm --exclude p3::max_total_ni --tag combined_v1
  # reproduce the original v2 files (19 SHOC+P3 parameters, plain lists)
  python make_json_files_multires_append.py --schemes shoc,p3 --tag v2 --legacy-format --force

Output files carry their parameter names, which create_ensembles.py reads:
  {"parameters": ["zm::tau", ...], "samples": [[...], ...]}
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

import ppe_parameters as reg

N_NE256 = 128
N_NE128 = 256
N_NE32 = 1024

LAMBDA_LOW, LAMBDA_HIGH = "shoc::lambda_low", "shoc::lambda_high"


def scale_value(u, p):
    if p.kind == "int":
        # equal-probability integer levels minimum..maximum
        lo, hi = int(p.minimum), int(p.maximum)
        return lo + min(int(u * (hi - lo + 1)), hi - lo)

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


def transform_samples(unit, params):
    keys = [p.key for p in params]
    swap = LAMBDA_LOW in keys and LAMBDA_HIGH in keys
    if swap:
        lo_idx, hi_idx = keys.index(LAMBDA_LOW), keys.index(LAMBDA_HIGH)

    samples = []
    for row in unit:
        vals = [scale_value(float(u), p) for u, p in zip(row, params)]

        # keep lambda_low <= lambda_high when both are sampled
        if swap and vals[lo_idx] > vals[hi_idx]:
            vals[lo_idx], vals[hi_idx] = vals[hi_idx], vals[lo_idx]

        samples.append(vals)

    return samples, swap


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


def write_json(path, obj, force):
    path = Path(path)
    if path.exists() and not force:
        raise FileExistsError(f"{path} exists; pass --force to overwrite")
    with open(path, "w") as f:
        json.dump(obj, f)


def split_list(values):
    out = []
    for v in values or []:
        out += [x.strip() for x in v.split(",") if x.strip()]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--schemes", required=True,
                    help=f"Comma-separated schemes to sample, from {reg.SCHEMES}")
    ap.add_argument("--exclude", action="append",
                    help="Parameter key(s) to leave out (repeatable or comma-separated), e.g. p3::max_total_ni")
    ap.add_argument("--tag", required=True, help="Suffix for output file names, e.g. zm_v1")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--legacy-format", action="store_true",
                    help="Write plain list-of-lists (only for the original 19 SHOC+P3 parameters)")
    ap.add_argument("--force", action="store_true", help="Overwrite existing output files")
    args = ap.parse_args()

    params = reg.select(schemes=split_list([args.schemes]), exclude=split_list(args.exclude))
    if not params:
        raise SystemExit("No parameters selected")
    keys = [p.key for p in params]
    if args.legacy_format and keys != reg.LEGACY_KEYS:
        raise SystemExit("--legacy-format needs exactly the original 19 SHOC+P3 parameters (--schemes shoc,p3, no --exclude)")

    out = Path(args.outdir)
    out_ne256 = out / f"append_nested_ne256_128_lhs_{args.tag}.json"
    out_ne128 = out / f"append_nested_ne128_256_lhs_{args.tag}.json"
    out_ne32 = out / f"append_nested_ne32_1024_lhs_{args.tag}.json"
    out_meta = out / f"append_nested_lhs_metadata_{args.tag}.json"

    ranks = build_appendable_ranks(
        n_full=N_NE32,
        n_mid=N_NE128,
        n_core=N_NE256,
        n_dim=len(params),
        seed=args.seed,
    )

    print("Rank-space nested LHS validation:")
    checks = validate_ranks(ranks)
    for n_subset, ok in checks:
        print(f"  first {n_subset:4d} rows are LHS: {ok}")
    if not all(ok for _, ok in checks):
        raise SystemExit("Nested LHS validation failed")

    unit = ranks_to_unit(ranks, args.seed)
    samples, swapped_lambda = transform_samples(unit, params)

    def payload(rows):
        return rows if args.legacy_format else {"parameters": keys, "samples": rows}

    write_json(out_ne256, payload(samples[:N_NE256]), args.force)
    write_json(out_ne128, payload(samples[:N_NE128]), args.force)
    write_json(out_ne32, payload(samples[:N_NE32]), args.force)

    metadata = {
        "method": "appendable_nested_lhs",
        "seed": args.seed,
        "tag": args.tag,
        "format": "list_of_lists" if args.legacy_format else "named",
        "n_ne256": N_NE256,
        "n_ne128": N_NE128,
        "n_ne32": N_NE32,
        "files": {
            "ne256": out_ne256.name,
            "ne128": out_ne128.name,
            "ne32": out_ne32.name,
        },
        "nesting": {
            "ne256_rows_in_ne32": [0, N_NE256 - 1],
            "ne128_rows_in_ne32": [0, N_NE128 - 1],
            "ne256_is_subset_of_ne128": True,
            "ne128_is_subset_of_ne32": True,
        },
        "schemes": split_list([args.schemes]),
        "excluded": split_list(args.exclude),
        "parameter_order": keys,
        "parameter_scales": [p.scale for p in params],
        "parameter_kinds": [p.kind for p in params],
        "parameter_ranges": [[p.minimum, p.maximum] for p in params],
        "parameter_defaults": [p.default for p in params],
        "enforce_lambda_low_le_lambda_high": swapped_lambda,
        "note": (
            "Strict LHS is validated in rank space before optional lambda_low/lambda_high swapping. "
            "If lambda ordering is enforced, lambda_low and lambda_high marginals may no longer be strict LHS. "
            "Integer parameters use equal-probability levels between minimum and maximum."
        ),
    }

    write_json(out_meta, metadata, args.force)

    print("\nWrote:")
    for f in (out_ne256, out_ne128, out_ne32, out_meta):
        print(f"  {f}")

    print("\nAppend structure:")
    print("  rows 0-127   -> ne256")
    print("  rows 0-255   -> ne128")
    print("  rows 0-1023  -> ne32")

    print("\nParameter order:")
    for i, p in enumerate(params, start=1):
        print(f"{i:2d}. {p.key:42s} {p.kind:4s} {p.scale:6s} min={p.minimum:g} max={p.maximum:g}")


if __name__ == "__main__":
    main()
