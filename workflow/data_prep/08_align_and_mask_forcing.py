# -*- coding: utf-8 -*-
import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import rasterio

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import (
    basin_paths,
    build_workspace_paths,
    default_config_path,
    ensure_workspace_dirs,
    legacy_file,
    load_legacy_module,
    read_config,
)


DATE_RE = re.compile(r"(\d{4}\.\d{2}\.\d{2})")


def collect_files(input_dir, start_date, end_date):
    records = []
    for file in Path(input_dir).glob("*.tif"):
        match = DATE_RE.search(file.name)
        if not match:
            continue
        date = datetime.strptime(match.group(1), "%Y.%m.%d")
        if start_date <= date <= end_date:
            records.append((date, file))
    records.sort(key=lambda item: item[0])
    return records


def write_series(module, input_dir, output_dir, prefix, dem_profile, dem_shape, basin_shape, start_date, end_date):
    records = collect_files(input_dir, start_date, end_date)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    count = 0
    for index, (date, input_file) in enumerate(records):
        output_file = Path(output_dir) / f"{index}_{prefix}_{date.strftime('%Y.%m.%d')}.tif"
        if output_file.exists():
            count += 1
            continue
        data = module.align_raster_to_dem(str(input_file), dem_profile, dem_shape)
        data = module.apply_basin_mask(data, basin_shape, dem_profile)
        profile = dem_profile.copy()
        profile.update(dtype=rasterio.float32, count=1, nodata=-9999, compress="lzw")
        with rasterio.open(output_file, "w", **profile) as dst:
            dst.write(data.astype("float32"), 1)
        count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description="Align precipitation, temperature, and ET rasters to the DEM grid.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--prec-source", choices=["mswep", "cmfd"], default="mswep")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    basin_info = basin_paths(config)
    ensure_workspace_dirs(paths)

    module = load_legacy_module(legacy_file("preprocess", "align_and_mask_rasters.py"))

    time_cfg = config["time"]
    start_date = datetime.strptime(f"{int(time_cfg['start_year'])}-01-01", "%Y-%m-%d")
    end_date = datetime.strptime(f"{int(time_cfg['end_year'])}-12-31", "%Y-%m-%d")

    with rasterio.open(paths["gis_dir"] / "dem_1km.tif") as dem:
        dem_profile = dem.profile.copy()
        dem_profile.update(dtype=rasterio.float32, nodata=-9999)
        dem_shape = (dem.height, dem.width)

    basin = gpd.read_file(basin_info["basin_shp"])
    basin_shape = [basin.geometry.iloc[0]]

    prec_input = paths["raw_prec_daily_dir"] if args.prec_source == "mswep" else paths["raw_prec_cmfd_daily_dir"]
    prec_output = paths["aligned_prec_dir"] if args.prec_source == "mswep" else paths["aligned_prec_cmfd_dir"]

    summary = {
        "precipitation": write_series(
            module,
            prec_input,
            prec_output,
            "PREC",
            dem_profile,
            dem_shape,
            basin_shape,
            start_date,
            end_date,
        ),
        "temperature": write_series(
            module,
            paths["raw_temp_daily_dir"],
            paths["aligned_temp_dir"],
            "TEMP",
            dem_profile,
            dem_shape,
            basin_shape,
            start_date,
            end_date,
        ),
        "evaporation": write_series(
            module,
            paths["raw_evap_daily_dir"],
            paths["aligned_evap_dir"],
            "EVAP",
            dem_profile,
            dem_shape,
            basin_shape,
            start_date,
            end_date,
        ),
    }

    for name, count in summary.items():
        print(f"{name}: {count} files")


if __name__ == "__main__":
    main()
