#!/usr/bin/env python3
"""
Create one PPE member case per JSON row by running a template script.

JSON formats:
  named   {"parameters": ["zm::tau", ...], "samples": [[...], ...]}   (from make_json_files_multires_append.py)
  legacy  [[...19 values...], ...]  columns in ppe_parameters.LEGACY_KEYS order (SHOC + P3)

Each member's template run gets:
  MEMBER_ID      m000, m001, ...
  PPE_CHANGES    "key=value key=value ..." for every varied parameter (run.ZM... applies these)
  SHOC_*/P3_*    the old per-parameter variables, for legacy keys only (old run.forENS... templates)
  do_case_build=false and OLD_EXECUTABLE when --old-exe is given

Parameters that are excluded, or not in the JSON, are not changed and run at their code default.

Examples:
  RES=ne32 python create_ensembles.py --json ./append_nested_ne32_1024_lhs_zm_v1.json \
      --template ./run.ZM.082019.F2010-SCREAMv1.pm-gpu.sh --old-exe <control>/build --start 0 --count 32
  ... --exclude-scheme p3 --exclude zm::mcsp_q_coeff     # hold these at default
  ... --dry-run                                         # validate and print changes, create nothing
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import ppe_parameters as reg

# Kept for tools that import the original 19-column order from here
PARAMS = [(p.legacy_env, p.key) for p in reg.PARAMETERS if p.legacy_env]


def fmt(v):
    # stable numeric formatting, preserves scientific notation when needed
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return format(v, ".17g")
    return str(v)


def run_one(template_path: Path, env: dict, log_path: Path, dry_run: bool):
    cmd = ["bash", str(template_path)]
    if dry_run:
        print(f"DRY RUN {env['MEMBER_ID']}: {env['PPE_CHANGES']}")
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("wb") as f:
        p = subprocess.Popen(cmd, env=env, stdout=f, stderr=subprocess.STDOUT)
        rc = p.wait()
    if rc != 0:
        raise RuntimeError(f"Member {env.get('MEMBER_ID')} failed rc={rc}, see {log_path}")


def split_list(values):
    out = []
    for v in values or []:
        out += [x.strip() for x in v.split(",") if x.strip()]
    return out


def load_json(path):
    """Return (format, parameter keys, rows)."""
    data = json.loads(Path(path).read_text())
    if isinstance(data, dict):
        if "parameters" not in data or "samples" not in data:
            raise ValueError('Named JSON must have "parameters" and "samples"')
        keys, rows = data["parameters"], data["samples"]
        fmt_name = "named"
    elif isinstance(data, list):
        keys, rows = list(reg.LEGACY_KEYS), data
        fmt_name = "legacy"
    else:
        raise ValueError("Expected a JSON list of rows or an object with parameters/samples")

    unknown = [k for k in keys if k not in reg.BY_KEY]
    if unknown:
        raise ValueError(f"Parameters not in ppe_parameters.py: {unknown}")
    dupes = sorted({k for k in keys if keys.count(k) > 1})
    if dupes:
        raise ValueError(f"Duplicate parameters in JSON: {dupes}")
    for i, row in enumerate(rows):
        if not isinstance(row, list) or len(row) != len(keys):
            raise ValueError(f"Row {i} must be a list of length {len(keys)}, "
                             f"got {type(row).__name__} len={len(row) if isinstance(row, list) else 'n/a'}")
        for k, v in zip(keys, row):
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                raise ValueError(f"Row {i} {k}: value {v!r} is not a number")
    return fmt_name, keys, rows


def coerce(p, v, idx):
    if p.kind == "int":
        if float(v) != int(v):
            raise ValueError(f"Row {idx} {p.key}: {v} is not an integer")
        return int(v)
    return float(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", required=True, help="Named or legacy (list-of-19) JSON file")
    ap.add_argument("--template", required=True, help="Path to run_template.sh")
    ap.add_argument("--old-exe", default="", help="Build directory of the control case to reuse (recommended)")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--count", type=int, default=0, help="0 means all")
    ap.add_argument("--exclude", action="append",
                    help="Parameter key(s) to hold at default (repeatable or comma-separated)")
    ap.add_argument("--exclude-scheme", action="append",
                    help=f"Scheme(s) to hold at default (repeatable or comma-separated), from {reg.SCHEMES}")
    ap.add_argument("--allow-out-of-range", action="store_true",
                    help="Warn instead of failing when a named JSON has values outside ppe_parameters.py ranges")
    ap.add_argument("--dry-run", action="store_true", help="Validate and print each member's changes, create nothing")
    args = ap.parse_args()

    fmt_name, keys, data = load_json(args.json)

    # Exclusions
    excl_keys = set(split_list(args.exclude))
    excl_schemes = set(split_list(args.exclude_scheme))
    bad = (excl_keys - set(reg.BY_KEY)) | {f"{s} (scheme)" for s in excl_schemes - set(reg.SCHEMES)}
    if bad:
        raise ValueError(f"Unknown exclusion(s): {sorted(bad)}")
    for k in sorted(excl_keys - set(keys)):
        print(f"WARNING: --exclude {k} is not in {args.json}; nothing to exclude")
    held = [k for k in keys if k in excl_keys or reg.BY_KEY[k].scheme in excl_schemes]
    active = [k for k in keys if k not in held]
    if not active:
        raise ValueError("Every parameter in the JSON is excluded; members would all be copies of the control")
    cols = [keys.index(k) for k in active]

    rows = data[args.start:] if args.count == 0 else data[args.start:args.start + args.count]
    if not rows:
        raise ValueError(f"No rows in range --start {args.start} --count {args.count} (file has {len(data)})")

    # Values: coerce types and check ranges for the whole selected range before creating anything
    members = []
    out_of_range = {}
    for idx, row in enumerate(rows, start=args.start):
        vals = []
        for k, c in zip(active, cols):
            p = reg.BY_KEY[k]
            v = coerce(p, row[c], idx)
            if not (p.minimum <= v <= p.maximum):
                out_of_range.setdefault(k, []).append((idx, v))
            vals.append(v)
        members.append((idx, vals))
    if out_of_range:
        lines = [f"  {k}: {len(v)} value(s) outside [{reg.BY_KEY[k].minimum:g}, {reg.BY_KEY[k].maximum:g}], "
                 f"e.g. row {v[0][0]} = {v[0][1]:g}" for k, v in out_of_range.items()]
        msg = "Values outside the ranges in ppe_parameters.py:\n" + "\n".join(lines)
        # legacy files predate the registry and used other ranges in earlier designs
        if fmt_name == "named" and not args.allow_out_of_range:
            raise ValueError(msg + "\nPass --allow-out-of-range if this is intended.")
        print("WARNING: " + msg)

    template = Path(args.template).resolve()
    if not template.is_file():
        raise FileNotFoundError(template)

    old_exe = Path(args.old_exe).resolve() if args.old_exe else None

    print(f"JSON: {args.json} ({fmt_name}, {len(data)} rows)")
    print(f"Members: m{args.start:03d}..m{members[-1][0]:03d} ({len(members)})")
    print(f"Varying {len(active)}: {' '.join(active)}")
    if held:
        print(f"Held at default {len(held)}: {' '.join(held)}")
    sys.stdout.flush()

    legacy_env_keys = [p.legacy_env for p in reg.PARAMETERS if p.legacy_env]
    for idx, vals in members:
        member_id = f"m{idx:03d}"
        env = os.environ.copy()
        # never inherit stale parameter settings from the calling shell
        for k in legacy_env_keys + ["PPE_CHANGES"]:
            env.pop(k, None)
        env["MEMBER_ID"] = member_id

        # Reuse executable if provided, skip build
        if old_exe:
            env["do_case_build"] = "false"
            env["OLD_EXECUTABLE"] = str(old_exe)

        env["PPE_CHANGES"] = " ".join(f"{k}={fmt(v)}" for k, v in zip(active, vals))
        for k, v in zip(active, vals):
            if reg.BY_KEY[k].legacy_env:
                env[reg.BY_KEY[k].legacy_env] = fmt(v)

        log_path = Path("ensemble_logs") / f"{member_id}.log"
        if not args.dry_run:
            print(f"Submitting {member_id}, log={log_path}")
            sys.stdout.flush()
        run_one(template, env, log_path, args.dry_run)


if __name__ == "__main__":
    main()
