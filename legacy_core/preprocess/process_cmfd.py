# -*- coding: utf-8 -*-
"""
Process CMFD V2.0 daily precipitation NetCDF into GeoTIFFs.

Input:
  .\data\raw\precipitation\CMFD\prec_CMFD_V0200_B-01_01dy_010deg_YYYYMM-YYYYMM.nc
Output:
  .\data\raw\precipitation\cmfd_daily\P_CMFD_YYYY.MM.DD.tif
"""
import os
import glob
import argparse
from datetime import datetime

import numpy as np
import xarray as xr
import rasterio
from rasterio.transform import from_origin

PROJECT_ROOT = r"."
CMFD_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "precipitation", "CMFD")
OUT_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "precipitation", "cmfd_daily")
FLOW_ACC_PATH = os.path.join(PROJECT_ROOT, "data", "gis", "flow_accumulation_masked.tif")

NODATA = -9999.0


def parse_args():
    parser = argparse.ArgumentParser(description="Process CMFD V2.0 precipitation to daily GeoTIFFs.")
    parser.add_argument("--start-date", type=str, default=None, help="Start date YYYY-MM-DD")
    parser.add_argument("--end-date", type=str, default=None, help="End date YYYY-MM-DD")
    parser.add_argument("--overwrite", action="store_true", default=False)
    return parser.parse_args()


def load_basin_bounds():
    with rasterio.open(FLOW_ACC_PATH) as src:
        bounds = src.bounds
    return bounds.left, bounds.bottom, bounds.right, bounds.top


def build_transform(lats, lons):
    lat_res = abs(float(lats[1] - lats[0])) if len(lats) > 1 else 0.1
    lon_res = abs(float(lons[1] - lons[0])) if len(lons) > 1 else 0.1
    west = float(lons.min()) - lon_res / 2.0
    north = float(lats.max()) + lat_res / 2.0
    return from_origin(west, north, lon_res, lat_res)


def main():
    args = parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    lon_min, lat_min, lon_max, lat_max = load_basin_bounds()
    start_date = datetime.strptime(args.start_date, "%Y-%m-%d").date() if args.start_date else None
    end_date = datetime.strptime(args.end_date, "%Y-%m-%d").date() if args.end_date else None

    nc_files = sorted(glob.glob(os.path.join(CMFD_DIR, "prec_CMFD_V0200_B-01_01dy_010deg_*.nc")))
    if not nc_files:
        print(f"[ERROR] No CMFD files found in: {CMFD_DIR}")
        return

    for nc_path in nc_files:
        ds = xr.open_dataset(nc_path)
        if "prec" not in ds.variables:
            ds.close()
            print(f"[WARN] Missing 'prec' in {os.path.basename(nc_path)}")
            continue

        prec = ds["prec"]
        prec_sub = prec.sel(
            lon=slice(lon_min, lon_max),
            lat=slice(lat_min, lat_max)
        )

        lons = prec_sub["lon"].values
        lats = prec_sub["lat"].values

        # Ensure north-up orientation (lat descending)
        if lats[0] < lats[-1]:
            lats = lats[::-1]
            prec_sub = prec_sub[:, ::-1, :]

        transform = build_transform(lats, lons)

        for i, t in enumerate(prec_sub["time"].values):
            date = datetime.utcfromtimestamp(np.datetime64(t, "s").astype(int)).date()
            if start_date and date < start_date:
                continue
            if end_date and date > end_date:
                continue

            out_name = f"P_CMFD_{date.strftime('%Y.%m.%d')}.tif"
            out_path = os.path.join(OUT_DIR, out_name)
            if (not args.overwrite) and os.path.exists(out_path):
                continue

            data = prec_sub[i].values.astype(np.float32) * 86400.0  # kg m-2 s-1 -> mm/day
            data = np.where(np.isfinite(data), data, NODATA)
            data = np.where(data < 0, 0.0, data)

            profile = {
                "driver": "GTiff",
                "height": data.shape[0],
                "width": data.shape[1],
                "count": 1,
                "dtype": "float32",
                "crs": "EPSG:4326",
                "transform": transform,
                "nodata": NODATA,
                "compress": "lzw",
            }

            with rasterio.open(out_path, "w", **profile) as dst:
                dst.write(data, 1)

        ds.close()
        print(f"[OK] {os.path.basename(nc_path)}")


if __name__ == "__main__":
    main()

