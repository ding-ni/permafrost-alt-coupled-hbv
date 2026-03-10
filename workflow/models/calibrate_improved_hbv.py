# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import (
    basin_paths,
    build_workspace_paths,
    default_config_path,
    ensure_workspace_dirs,
    legacy_file,
    load_legacy_module,
    observation_mode_for_legacy,
    patch_module,
    read_config,
    temporary_argv,
)


def main():
    parser = argparse.ArgumentParser(description="Calibrate the improved distributed HBV model.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--prec-source", choices=["mswep", "cmfd"], default=None)
    parser.add_argument("--glacier-mode", choices=["inline", "series", "off"], default="inline")
    parser.add_argument("--maxiter", type=int, default=80)
    parser.add_argument("--popsize", type=int, default=12)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--quick-test", action="store_true")
    parser.add_argument("--quick-days", type=int, default=30)
    parser.add_argument("--skip-gm", action="store_true")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    basin = basin_paths(config)
    ensure_workspace_dirs(paths)
    time_cfg = config["time"]
    prec_source = args.prec_source or config.get("default_precip_source", "mswep")

    module = load_legacy_module(legacy_file("models", "improved_hbv.py"))
    patch_module(
        module,
        {
            "PROJECT_ROOT": str(paths["workspace_root"]),
            "PREC_DIR_MSWEP": str(paths["aligned_prec_dir"]),
            "PREC_DIR_CMFD": str(paths["aligned_prec_cmfd_dir"]),
            "PREC_DIR": str(paths["aligned_prec_dir"] if prec_source == "mswep" else paths["aligned_prec_cmfd_dir"]),
            "TEMP_DIR": str(paths["aligned_temp_dir"]),
            "EVAP_DIR": str(paths["aligned_evap_dir"]),
            "OBS_FILE": str(basin["obs_csv"]),
            "FLOW_ACC_PATH": str(paths["gis_dir"] / "flow_accumulation_masked.tif"),
            "GLACIER_MELT_DIR": str(paths["glacier_melt_dir"]),
            "GLACIER_MASK_PATH": str(paths["gis_dir"] / "glacier_mask.tif"),
            "LOG_DIR": str(paths["logs_dir"]),
            "CACHE_DIR": str(paths["cache_dir"]),
            "CATCHMENT_AREA": float(config["catchment_area_km2"]),
            "WARMUP_START": time_cfg["warmup_start"],
            "WARMUP_END": time_cfg["warmup_end"],
            "CALIB_START": time_cfg["calibration_start"],
            "CALIB_END": time_cfg["calibration_end"],
            "VALID_START": time_cfg["validation_start"],
            "VALID_END": time_cfg["validation_end"],
            "SIM_START": time_cfg["calibration_start"],
            "SIM_END": time_cfg["validation_end"],
            "OBS_MODE_OVERRIDE": observation_mode_for_legacy(config),
        },
    )

    argv = [
        "calibrate_improved_hbv.py",
        "--maxiter",
        str(args.maxiter),
        "--popsize",
        str(args.popsize),
        "--seed",
        str(args.seed),
        "--workers",
        str(args.workers),
        "--prec-source",
        prec_source,
        "--glacier-mode",
        args.glacier_mode,
    ]
    if args.quick_test:
        argv.extend(["--quick-test", "--quick-days", str(args.quick_days)])
    if args.skip_gm:
        argv.append("--skip-gm")

    with temporary_argv(argv):
        module.main()


if __name__ == "__main__":
    main()
