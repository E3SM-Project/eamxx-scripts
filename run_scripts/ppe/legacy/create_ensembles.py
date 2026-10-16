#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
from pathlib import Path

PARAMS = [
    ("SHOC_THL2TUNE", "shoc::thl2tune"),
    ("SHOC_QW2TUNE", "shoc::qw2tune"),
    ("SHOC_LENGTH_FAC", "shoc::length_fac"),
    ("SHOC_C_DIAG_3RD_MOM", "shoc::c_diag_3rd_mom"),
    ("SHOC_COEFF_KH", "shoc::coeff_kh"),
    ("SHOC_COEFF_KM", "shoc::coeff_km"),
    ("SHOC_LAMBDA_LOW", "shoc::lambda_low"),
    ("SHOC_LAMBDA_HIGH", "shoc::lambda_high"),
    ("P3_SPA_CCN_TO_NC_FACTOR", "p3::spa_ccn_to_nc_factor"),
    ("P3_CLDLIQ_TO_ICE_COLLECTION_FACTOR", "p3::cldliq_to_ice_collection_factor"),
    ("P3_RAIN_TO_ICE_COLLECTION_FACTOR", "p3::rain_to_ice_collection_factor"),
    ("P3_ACCRETION_PREFACTOR", "p3::accretion_prefactor"),
    ("P3_DEPOSITION_NUCLEATION_EXPONENT", "p3::deposition_nucleation_exponent"),
    ("P3_MAX_TOTAL_NI", "p3::max_total_ni"),
    ("P3_ICE_SEDIMENTATION_FACTOR", "p3::ice_sedimentation_factor"),
    ("P3_RAIN_SELFCOLLECTION_BREAKUP_DIAMETER", "p3::rain_selfcollection_breakup_diameter"),
    ("P3_AUTOCONVERSION_PREFACTOR", "p3::autoconversion_prefactor"),
    ("P3_AUTOCONVERSION_QC_EXPONENT", "p3::autoconversion_qc_exponent"),
    ("P3_AUTOCONVERSION_RADIUS", "p3::autoconversion_radius"),
]

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
        print("DRY RUN:", " ".join(cmd), "MEMBER_ID=", env.get("MEMBER_ID", ""))
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("wb") as f:
        p = subprocess.Popen(cmd, env=env, stdout=f, stderr=subprocess.STDOUT)
        rc = p.wait()
    if rc != 0:
        raise RuntimeError(f"Member {env.get('MEMBER_ID')} failed rc={rc}, see {log_path}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True, help="Path to JSON file containing [[...19...],[...],...]")
    ap.add_argument("--template", required=True, help="Path to run_template.sh")
    ap.add_argument("--old-exe", default="", help="Path to prebuilt e3sm.exe to reuse (recommended)")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--count", type=int, default=0, help="0 means all")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    data = json.loads(Path(args.json).read_text())
    if not isinstance(data, list):
        raise ValueError("Expected top-level JSON array of rows (list-of-lists).")

    rows = data[args.start:] if args.count == 0 else data[args.start:args.start + args.count]

    template = Path(args.template).resolve()
    if not template.is_file():
        raise FileNotFoundError(template)

    old_exe = Path(args.old_exe).resolve() if args.old_exe else None
    #if old_exe and not old_exe.is_file():
    #    raise FileNotFoundError(old_exe)

    for idx, row in enumerate(rows, start=args.start):
        if not isinstance(row, list) or len(row) != len(PARAMS):
            raise ValueError(f"Row {idx} must be a list of length {len(PARAMS)}, got {type(row)} len={len(row) if isinstance(row, list) else 'n/a'}")

        member_id = f"m{idx:03d}"
        env = os.environ.copy()
        env["MEMBER_ID"] = member_id

        # Reuse executable if provided, skip build
        if old_exe:
            env["do_case_build"] = "false"
            env["OLD_EXECUTABLE"] = str(old_exe)

        # Populate the 19 tuning env vars
        for (env_key, _atm_key), value in zip(PARAMS, row):
            env[env_key] = fmt(value)

        log_path = Path("ensemble_logs") / f"{member_id}.log"
        print(f"Submitting {member_id}, log={log_path}")
        run_one(template, env, log_path, args.dry_run)

if __name__ == "__main__":
    main()
