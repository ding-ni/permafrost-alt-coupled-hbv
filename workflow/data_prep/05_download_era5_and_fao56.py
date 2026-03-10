# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import (
    bbox_list,
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
    parser = argparse.ArgumentParser(description="Download ERA5-Land and FAO-56 support variables.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--content", choices=["era5", "fao56", "all", "mswep-notes"], default="all")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    ensure_workspace_dirs(paths)
    years = list(year_range(config))

    if args.content in {"era5", "all", "mswep-notes"}:
        module = load_legacy_module(legacy_file("preprocess", "download_meteorological_data.py"))
        patch_module(
            module,
            {
                "PROJECT_ROOT": str(paths["workspace_root"]),
                "RAW_TEMP_DIR": str(paths["raw_temp_dir"]),
                "RAW_EVAP_DIR": str(paths["raw_evap_dir"]),
                "RAW_PREC_DIR": str(paths["raw_prec_root"]),
                "BASIN_BBOX": bbox_list(config),
                "TUOTUOHE_BBOX": bbox_list(config),
                "START_YEAR": years[0],
                "END_YEAR": years[-1],
            },
        )
        if args.content in {"era5", "all"}:
            module.download_era5_all()
        if args.content in {"mswep-notes", "all"}:
            module.print_mswep_instructions()

    if args.content in {"fao56", "all"}:
        module = load_legacy_module(legacy_file("preprocess", "download_fao56_variables.py"))
        patch_module(
            module,
            {
                "PROJECT_ROOT": str(paths["workspace_root"]),
                "BASIN_BBOX": bbox_list(config),
                "TUOTUOHE_BBOX": bbox_list(config),
                "START_YEAR": years[0],
                "END_YEAR": years[-1],
                "RAW_DIR": str(paths["raw_root"]),
                "SOLAR_DIR": str(paths["raw_solar_dir"]),
                "WIND_DIR": str(paths["raw_wind_dir"]),
                "DEWPOINT_DIR": str(paths["raw_dewpoint_dir"]),
            },
        )
        module.main()


if __name__ == "__main__":
    main()
