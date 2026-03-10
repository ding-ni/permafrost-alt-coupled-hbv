# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import (
    build_workspace_paths,
    default_config_path,
    ensure_workspace_dirs,
    legacy_file,
    load_legacy_module,
    patch_module,
    read_config,
    year_range,
)


def main():
    parser = argparse.ArgumentParser(description="Process ERA5 temperature and evapotranspiration data.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--et-method", choices=["fao56", "era5"], default="fao56")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    ensure_workspace_dirs(paths)
    years = list(year_range(config))

    era5_module = load_legacy_module(legacy_file("preprocess", "process_era5_netcdf.py"))
    patch_module(
        era5_module,
        {
            "PROJECT_ROOT": str(paths["workspace_root"]),
            "RAW_TEMP_DIR": str(paths["raw_temp_dir"]),
            "RAW_EVAP_DIR": str(paths["raw_evap_dir"]),
            "TEMP_DAILY_DIR": str(paths["raw_temp_daily_dir"]),
            "EVAP_DAILY_DIR": str(paths["raw_evap_daily_dir"]),
            "YEARS": years,
        },
    )

    for year in years:
        era5_module.process_temperature(year)
        if args.et_method == "era5":
            era5_module.process_evaporation(year)

    if args.et_method == "fao56":
        et_module = load_legacy_module(legacy_file("preprocess", "calculate_fao56_et.py"))
        patch_module(
            et_module,
            {
                "PROJECT_ROOT": str(paths["workspace_root"]),
                "TEMP_DIR": str(paths["raw_temp_dir"]),
                "SOLAR_DIR": str(paths["raw_solar_dir"]),
                "WIND_DIR": str(paths["raw_wind_dir"]),
                "DEWPOINT_DIR": str(paths["raw_dewpoint_dir"]),
                "OUTPUT_DIR": str(paths["raw_evap_daily_dir"]),
                "YEARS": years,
                "ELEVATION": float(config.get("fao56_mean_elevation_m", 4500.0)),
            },
        )
        for year in years:
            et_module.process_year(year)


if __name__ == "__main__":
    main()
