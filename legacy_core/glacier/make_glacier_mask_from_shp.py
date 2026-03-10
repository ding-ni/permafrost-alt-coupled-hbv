# -*- coding: utf-8 -*-
"""
Create glacier mask raster from glacier shapefile using DEM grid.

Outputs:
  - data/gis/glacier_mask.tif
"""
import os
import sys

sys.path.insert(0, r".\Hapi-main\Hapi-main\src")

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from math import sin, radians

PROJECT_ROOT = r"."
DEM_FILE = os.path.join(PROJECT_ROOT, "data", "gis", "dem_1km.tif")
BASIN_SHP = os.path.join(PROJECT_ROOT, "data", "gis", "basin.shp")
GLACIER_SHP = os.path.join(PROJECT_ROOT, "data", "gis", "glacier_shp", "glacier.shp")
OUT_FILE = os.path.join(PROJECT_ROOT, "data", "gis", "glacier_mask.tif")


def compute_area_km2_from_mask(mask, transform):
    """Compute area for EPSG:4326 rasters."""
    pixel_width_deg = transform.a
    pixel_height_deg = -transform.e
    top_lat = transform.f
    r = 6371000.0
    dlon = radians(pixel_width_deg)
    total_area_m2 = 0.0
    for row in range(mask.shape[0]):
        count = int((mask[row, :] > 0).sum())
        if count == 0:
            continue
        lat_n = top_lat - row * pixel_height_deg
        lat_s = top_lat - (row + 1) * pixel_height_deg
        area_row_m2 = (r ** 2) * dlon * (sin(radians(lat_n)) - sin(radians(lat_s)))
        total_area_m2 += area_row_m2 * count
    return total_area_m2 / 1e6


def main():
    if not os.path.exists(DEM_FILE):
        raise FileNotFoundError(DEM_FILE)
    if not os.path.exists(GLACIER_SHP):
        raise FileNotFoundError(GLACIER_SHP)

    with rasterio.open(DEM_FILE) as dem:
        profile = dem.profile.copy()
        transform = dem.transform
        out_shape = (dem.height, dem.width)
        dem_crs = dem.crs

    glacier_gdf = gpd.read_file(GLACIER_SHP)
    if glacier_gdf.crs != dem_crs:
        glacier_gdf = glacier_gdf.to_crs(dem_crs)

    glacier_shapes = [
        (geom, 1) for geom in glacier_gdf.geometry if geom is not None and not geom.is_empty
    ]
    if not glacier_shapes:
        raise ValueError("No valid glacier geometries found.")

    glacier_mask = rasterize(
        glacier_shapes,
        out_shape=out_shape,
        transform=transform,
        fill=0,
        dtype="uint8",
    )

    if os.path.exists(BASIN_SHP):
        basin_gdf = gpd.read_file(BASIN_SHP)
        if basin_gdf.crs != dem_crs:
            basin_gdf = basin_gdf.to_crs(dem_crs)
        basin_shapes = [
            (geom, 1) for geom in basin_gdf.geometry if geom is not None and not geom.is_empty
        ]
        basin_mask = rasterize(
            basin_shapes,
            out_shape=out_shape,
            transform=transform,
            fill=0,
            dtype="uint8",
        )
        glacier_mask = np.where(basin_mask == 1, glacier_mask, 255).astype("uint8")
        nodata = 255
    else:
        nodata = 0

    profile.update(dtype="uint8", count=1, nodata=nodata, compress="lzw")
    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    with rasterio.open(OUT_FILE, "w", **profile) as dst:
        dst.write(glacier_mask, 1)

    valid = (glacier_mask > 0) & (glacier_mask != nodata)
    area_km2 = compute_area_km2_from_mask(valid.astype("uint8"), transform)
    print(f"[OK] glacier_mask.tif written: {OUT_FILE}")
    print(f"      glacier pixels: {int(valid.sum())}")
    print(f"      glacier area: {area_km2:.2f} km2")


if __name__ == "__main__":
    main()

