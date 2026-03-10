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
    parser = argparse.ArgumentParser(description="Calibrate and run the original distributed HBV baseline.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--maxiter", type=int, default=80)
    parser.add_argument("--popsize", type=int, default=12)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    basin = basin_paths(config)
    ensure_workspace_dirs(paths)
    time_cfg = config["time"]

    module = load_legacy_module(legacy_file("models", "original_hbv.py"))
    patch_module(
        module,
        {
            "PROJECT_ROOT": paths["workspace_root"],
            "THESIS_ROOT": paths["workspace_root"],
            "PREC_DIR": paths["aligned_prec_dir"],
            "TEMP_DIR": paths["aligned_temp_dir"],
            "EVAP_DIR": paths["aligned_evap_dir"],
            "OBS_FILE": basin["obs_csv"],
            "FLOW_ACC_PATH": paths["gis_dir"] / "flow_accumulation_masked.tif",
            "RESULTS_DIR": paths["original_results_dir"],
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
        "run_original_hbv.py",
        "--maxiter",
        str(args.maxiter),
        "--popsize",
        str(args.popsize),
        "--seed",
        str(args.seed),
    ]
    with temporary_argv(argv):
        module.main()


if __name__ == "__main__":
    main()
