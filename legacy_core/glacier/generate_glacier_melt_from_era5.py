# -*- coding: utf-8 -*-
"""
Generate daily glacier melt rasters (mm/d) using ERA5 ssrd + temperature.

Outputs:
  data/aligned_masked/glacier_melt/GM_YYYY.MM.DD.tif
"""
import argparse
import os
import re
from datetime import datetime

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
import xarray as xr


PROJECT_ROOT = r"."
TEMP_DIR = os.path.join(PROJECT_ROOT, "data", "aligned_masked", "temp")
RAD_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "solar_radiation")
OUT_DIR = os.path.join(PROJECT_ROOT, "data", "aligned_masked", "glacier_melt")
DEM_PATH = os.path.join(PROJECT_ROOT, "data", "gis", "dem_1km.tif")
FLOW_ACC_PATH = os.path.join(PROJECT_ROOT, "data", "gis", "flow_accumulation_masked.tif")
GLACIER_MASK_PATH = os.path.join(PROJECT_ROOT, "data", "gis", "glacier_mask.tif")

# Glacier melt parameters
TTM = 0.0           # 温度融化阈值 (°C)
TTM_RAD = -5.0      # 辐射融化阈值 (°C) - 低于此温度不考虑辐射融化
CFMAX_GLACIER = 5.5 # 度日因子 (mm/°C/d)
RAD_COEF = 0.002    # 辐射系数 (降低以减少辐射融化贡献)

NODATA = -9999.0
OVERWRITE = True


def parse_args():
    parser = argparse.ArgumentParser(description="Generate glacier melt rasters from ERA5 ssrd.")
    parser.add_argument("--start-year", type=int, default=None, help="Start year (inclusive).")
    parser.add_argument("--end-year", type=int, default=None, help="End year (inclusive).")
    parser.add_argument("--overwrite", dest="overwrite", action="store_true")
    parser.add_argument("--no-overwrite", dest="overwrite", action="store_false")
    parser.set_defaults(overwrite=OVERWRITE)
    return parser.parse_args()


def parse_date_from_name(name):
    match = re.search(r"\d{4}\.\d{2}\.\d{2}", name)
    if not match:
        return None
    return datetime.strptime(match.group(0), "%Y.%m.%d").date()


def build_src_transform(lats, lons):
    lat_res = abs(float(lats[1] - lats[0]))
    lon_res = abs(float(lons[1] - lons[0]))
    west = float(lons.min()) - lon_res / 2.0
    north = float(lats.max()) + lat_res / 2.0
    return from_origin(west, north, lon_res, lat_res)


def load_daily_radiation(year):
    nc_path = os.path.join(RAD_DIR, f"era5_ssrd_{year}.nc")
    if not os.path.exists(nc_path):
        return None, None, None

    ds = xr.open_dataset(nc_path)
    if "ssrd" not in ds.variables:
        ds.close()
        return None, None, None

    ssrd = ds["ssrd"]
    daily = ssrd.resample(valid_time="1D").sum()

    lats = ds["latitude"].values
    lons = ds["longitude"].values
    daily_dates = [pd.to_datetime(t).date() for t in daily["valid_time"].values]

    data = daily.values
    ds.close()

    return daily_dates, data, (lats, lons)


def main():
    args = parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    with rasterio.open(DEM_PATH) as dem:
        ref_profile = dem.profile.copy()
        ref_transform = dem.transform
        ref_crs = dem.crs
        ref_shape = dem.shape

    with rasterio.open(FLOW_ACC_PATH) as src:
        flow_acc = src.read(1)
        flow_nodata = src.nodata
    if flow_nodata is None:
        basin_mask = np.isfinite(flow_acc)
    else:
        basin_mask = flow_acc != flow_nodata

    with rasterio.open(GLACIER_MASK_PATH) as src:
        glacier_raw = src.read(1)
        glacier_nodata = src.nodata
    if glacier_nodata is not None:
        glacier_raw = np.where(glacier_raw == glacier_nodata, 0, glacier_raw)
    glacier_mask = glacier_raw > 0
    glacier_mask = glacier_mask & basin_mask

    temp_files = [f for f in os.listdir(TEMP_DIR) if f.lower().endswith(".tif")]
    temp_by_year = {}
    for fname in temp_files:
        date = parse_date_from_name(fname)
        if date is None:
            continue
        temp_by_year.setdefault(date.year, []).append((date, os.path.join(TEMP_DIR, fname)))

    years = sorted(temp_by_year.keys())
    if args.start_year is not None:
        years = [y for y in years if y >= args.start_year]
    if args.end_year is not None:
        years = [y for y in years if y <= args.end_year]
    print(f"Found temperature years: {years}")

    for year in years:
        dates, rad_data, (lats, lons) = load_daily_radiation(year)
        if dates is None:
            print(f"[WARN] Missing radiation for {year}, skipping")
            continue

        src_transform = build_src_transform(lats, lons)
        src_crs = "EPSG:4326"

        rad_by_date = {d: rad_data[i] for i, d in enumerate(dates)}

        year_files = sorted(temp_by_year[year], key=lambda x: x[0])
        print(f"{year}: {len(year_files)} days")

        for date, temp_path in year_files:
            out_name = f"GM_{date.strftime('%Y.%m.%d')}.tif"
            out_path = os.path.join(OUT_DIR, out_name)
            if (not args.overwrite) and os.path.exists(out_path):
                continue

            if date not in rad_by_date:
                continue

            with rasterio.open(temp_path) as src:
                temp = src.read(1)
                temp_nodata = src.nodata

            rad_src = rad_by_date[date].astype(np.float32)
            rad_dst = np.full(ref_shape, np.nan, dtype=np.float32)
            reproject(
                source=rad_src,
                destination=rad_dst,
                src_transform=src_transform,
                src_crs=src_crs,
                dst_transform=ref_transform,
                dst_crs=ref_crs,
                resampling=Resampling.bilinear,
            )

            if temp_nodata is not None:
                temp = np.where(temp == temp_nodata, np.nan, temp)

            rad_wm2 = np.where(np.isfinite(rad_dst), rad_dst, 0.0) / 86400.0

            melt = np.zeros(ref_shape, dtype=np.float32)
            temp_eff = np.where(np.isfinite(temp), temp, -999.0)

            # 温度融化: 仅当T > TTM时
            melt_temp = np.maximum(temp_eff - TTM, 0.0) * CFMAX_GLACIER

            # 辐射融化: 仅当T > TTM_RAD时 (避免冬季不合理融水)
            rad_melt = np.where(temp_eff > TTM_RAD,
                               np.maximum(rad_wm2, 0.0) * RAD_COEF,
                               0.0)

            melt_val = melt_temp + rad_melt
            melt_val = np.maximum(melt_val, 0.0)

            melt[glacier_mask] = melt_val[glacier_mask]
            melt[~basin_mask] = NODATA

            profile = ref_profile.copy()
            profile.update(dtype=rasterio.float32, count=1, nodata=NODATA)
            with rasterio.open(out_path, "w", **profile) as dst:
                dst.write(melt, 1)

        print(f"[OK] {year} done")


if __name__ == "__main__":
    main()

