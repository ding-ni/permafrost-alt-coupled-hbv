# -*- coding: utf-8 -*-
import argparse
import json
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.public_aliases import legacy_to_public_scheme_name, public_to_legacy_scheme_names
from shared.runtime import (
    basin_paths,
    build_workspace_paths,
    coupling_start_year,
    default_config_path,
    ensure_workspace_dirs,
    legacy_file,
    load_full_v6_params,
    load_legacy_module,
    observation_mode_for_legacy,
    patch_module,
    read_config,
    temporary_argv,
)


def build_legacy_v6_params(params):
    return {
        "TT": params["TT"],
        "FC": params["FC"],
        "BETA": params["BETA"],
        "LP": params["LP"],
        "RFCF": params["RFCF"],
        "SFCF": params["SFCF"],
        "CFR": params["CFR"],
        "CWH": params["CWH"],
        "CFMAX_mid": params["CFMAX_MID"],
        "CFMAX_high": params["CFMAX_HIGH"],
        "K": params["K"],
        "K0": params.get("K0", params["K"]),
        "K1": params["K1"],
        "K2": params["K2"],
        "PERC": params["PERC"],
        "K_MUSK": params["K_MUSK"],
        "ICE_FACTOR": params["ICE_FACTOR"],
    }


def rewrite_summary(summary_path):
    if not summary_path.exists():
        return
    with summary_path.open("r", encoding="utf-8") as handle:
        summary = json.load(handle)

    if "best_scheme" in summary:
        summary["best_scheme"] = legacy_to_public_scheme_name(summary["best_scheme"])
    if "best_scheme_obj" in summary:
        summary["best_scheme_obj"] = legacy_to_public_scheme_name(summary["best_scheme_obj"])
    if "best_scheme_full" in summary:
        summary["best_scheme_full"] = legacy_to_public_scheme_name(summary["best_scheme_full"])

    schemes = summary.get("schemes", {})
    renamed = {}
    for legacy_name, payload in schemes.items():
        public_name = legacy_to_public_scheme_name(legacy_name)
        renamed[public_name] = payload
    summary["schemes"] = renamed

    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)


def rename_scheme_outputs(run_dir):
    for file in run_dir.glob("*_timeseries.csv"):
        legacy_name = file.name[: -len("_timeseries.csv")]
        public_name = legacy_to_public_scheme_name(legacy_name)
        if public_name != legacy_name:
            shutil.move(str(file), str(run_dir / f"{public_name}_timeseries.csv"))
    for file in run_dir.glob("*_param_stats.csv"):
        legacy_name = file.name[: -len("_param_stats.csv")]
        public_name = legacy_to_public_scheme_name(legacy_name)
        if public_name != legacy_name:
            shutil.move(str(file), str(run_dir / f"{public_name}_param_stats.csv"))


def newest_run_dir(runs_dir, before):
    candidates = [path for path in runs_dir.iterdir() if path.is_dir() and path not in before]
    if not candidates:
        return None
    return max(candidates, key=lambda item: item.stat().st_mtime)


def main():
    parser = argparse.ArgumentParser(description="Run the ALT-coupled HBV multi-scheme comparison.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--improved-run", default="")
    parser.add_argument("--only", default="", help="Comma-separated public scheme names.")
    parser.add_argument("--alt-ref", default="", help="ALT reference raster name or absolute path.")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    basin = basin_paths(config)
    ensure_workspace_dirs(paths)
    time_cfg = config["time"]
    full_v6_params = build_legacy_v6_params(load_full_v6_params(config, args.improved_run or None))

    module = load_legacy_module(legacy_file("models", "alt_multi_scheme.py"))
    patch_module(
        module,
        {
            "PROJECT_ROOT": str(paths["workspace_root"]),
            "PREC_DIR": str(paths["aligned_prec_dir"]),
            "TEMP_DIR": str(paths["aligned_temp_dir"]),
            "EVAP_DIR": str(paths["aligned_evap_dir"]),
            "OBS_FILE": str(basin["obs_csv"]),
            "FLOW_ACC_PATH": str(paths["gis_dir"] / "flow_accumulation_masked.tif"),
            "GLACIER_MASK_PATH": str(paths["gis_dir"] / "glacier_mask.tif"),
            "ALT_DIR": str(paths["alt_basin_dir"]),
            "SF_DIR": str(paths["alt_scale_dir"]),
            "RESULTS_DIR": str(paths["runs_dir"]),
            "CATCHMENT_AREA": float(config["catchment_area_km2"]),
            "V6_PARAMS": full_v6_params,
            "WARMUP_START": time_cfg["warmup_start"],
            "WARMUP_END": time_cfg["warmup_end"],
            "CALIB_START": time_cfg["calibration_start"],
            "CALIB_END": time_cfg["calibration_end"],
            "VALID_START": time_cfg["validation_start"],
            "VALID_END": time_cfg["validation_end"],
            "SIM_START": time_cfg["calibration_start"],
            "SIM_END": time_cfg["validation_end"],
            "OBS_MODE_OVERRIDE": observation_mode_for_legacy(config),
            "COUPLING_START_YEAR": coupling_start_year(config),
        },
    )

    before = {path for path in paths["runs_dir"].glob("*") if path.is_dir()}
    argv = ["run_alt_coupled_hbv.py"]
    if args.only:
        legacy_names = public_to_legacy_scheme_names(args.only.split(","))
        argv.extend(["--only", ",".join(legacy_names)])
    if args.alt_ref:
        argv.extend(["--alt-ref", args.alt_ref])

    with temporary_argv(argv):
        module.main()

    run_dir = newest_run_dir(paths["runs_dir"], before)
    if run_dir is not None:
        rename_scheme_outputs(run_dir)
        rewrite_summary(run_dir / "summary.json")


if __name__ == "__main__":
    main()
