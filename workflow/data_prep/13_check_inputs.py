# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import basin_paths, build_workspace_paths, default_config_path, expected_days, read_config


def tif_count(path):
    path = Path(path)
    if not path.exists():
        return 0
    return len(list(path.glob("*.tif")))


def check_file(path):
    path = Path(path)
    exists = path.exists()
    print(f"{'OK' if exists else 'MISSING'}  {path}")
    return exists


def main():
    parser = argparse.ArgumentParser(description="Quickly check whether required basin inputs are present.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--prec-source", choices=["mswep", "cmfd"], default=None)
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    basin = basin_paths(config)
    prec_source = args.prec_source or config.get("default_precip_source", "mswep")
    expected = expected_days(config)

    print("=" * 72)
    print(f"Basin: {config['basin_name_en']}")
    print(f"Expected daily steps: {expected}")
    print(f"Precipitation source: {prec_source}")
    print("=" * 72)

    check_file(basin["basin_shp"])
    check_file(paths["gis_dir"] / "dem_1km.tif")
    check_file(paths["gis_dir"] / "flow_direction.tif")
    check_file(paths["gis_dir"] / "flow_accumulation.tif")
    check_file(paths["gis_dir"] / "flow_accumulation_masked.tif")
    check_file(paths["gis_dir"] / "elevation_zone_mid.tif")
    check_file(paths["gis_dir"] / "elevation_zone_high.tif")
    check_file(basin["obs_csv"])

    alt_count = tif_count(paths["alt_basin_dir"])
    global_param_count = tif_count(paths["global_param_dir"])
    alt_param_count = tif_count(paths["alt_param_dir"])
    temp_count = tif_count(paths["aligned_temp_dir"])
    evap_count = tif_count(paths["aligned_evap_dir"])
    prec_count = tif_count(paths["aligned_prec_dir"] if prec_source == "mswep" else paths["aligned_prec_cmfd_dir"])

    print("-" * 72)
    print(f"ALT rasters: {alt_count}")
    print(f"Global HBV parameter rasters: {global_param_count}")
    print(f"ALT-adjusted parameter rasters: {alt_param_count}")
    print(f"Precipitation rasters: {prec_count}")
    print(f"Temperature rasters: {temp_count}")
    print(f"Evapotranspiration rasters: {evap_count}")
    print("-" * 72)


if __name__ == "__main__":
    main()
