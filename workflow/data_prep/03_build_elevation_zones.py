# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

import numpy as np
import rasterio

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import build_workspace_paths, default_config_path, ensure_workspace_dirs, read_config


def save_zone(profile, data, output_file):
    profile = profile.copy()
    profile.update(dtype=rasterio.uint8, count=1, nodata=0, compress="lzw")
    with rasterio.open(output_file, "w", **profile) as dst:
        dst.write(data.astype(np.uint8), 1)


def main():
    parser = argparse.ArgumentParser(description="Build low/mid/high elevation zoning rasters.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--mid-elevation-min", type=float, default=4500.0)
    parser.add_argument("--high-elevation-min", type=float, default=5000.0)
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    ensure_workspace_dirs(paths)

    dem_file = paths["gis_dir"] / "dem_1km.tif"
    with rasterio.open(dem_file) as src:
        dem = src.read(1)
        profile = src.profile.copy()

    valid = np.isfinite(dem) & (dem > 0)
    zone_low = ((dem < args.mid_elevation_min) & valid).astype(np.uint8)
    zone_mid = ((dem >= args.mid_elevation_min) & (dem < args.high_elevation_min) & valid).astype(np.uint8)
    zone_high = ((dem >= args.high_elevation_min) & valid).astype(np.uint8)

    save_zone(profile, zone_low, paths["gis_dir"] / "elevation_zone_low.tif")
    save_zone(profile, zone_mid, paths["gis_dir"] / "elevation_zone_mid.tif")
    save_zone(profile, zone_high, paths["gis_dir"] / "elevation_zone_high.tif")
    save_zone(profile, zone_low + zone_mid * 2 + zone_high * 3, paths["gis_dir"] / "elevation_zones.tif")

    print(f"Low-elevation cells: {int(zone_low.sum())}")
    print(f"Mid-elevation cells: {int(zone_mid.sum())}")
    print(f"High-elevation cells: {int(zone_high.sum())}")
    print(f"Total valid cells: {int(zone_low.sum() + zone_mid.sum() + zone_high.sum())}")


if __name__ == "__main__":
    main()
