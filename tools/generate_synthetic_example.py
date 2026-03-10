# -*- coding: utf-8 -*-
import argparse
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Polygon


def write_raster(path, data, transform, crs="EPSG:4326", nodata=-9999.0):
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
        nodata=nodata,
        compress="lzw",
    ) as dst:
        dst.write(data.astype("float32"), 1)


def write_uint_raster(path, data, transform, crs="EPSG:4326", nodata=0):
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype="uint8",
        crs=crs,
        transform=transform,
        nodata=nodata,
        compress="lzw",
    ) as dst:
        dst.write(data.astype("uint8"), 1)


def main():
    parser = argparse.ArgumentParser(description="Generate a lightweight synthetic basin example.")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config_path = Path(args.config)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    repo_root = config_path.resolve().parents[2]

    workspace_root = (config_path.parent / config["workspace_root"]).resolve()
    inputs_root = config_path.parent / "inputs"
    alt_inputs_root = inputs_root / "alt"

    for path in [workspace_root, inputs_root, alt_inputs_root]:
        path.mkdir(parents=True, exist_ok=True)

    polygon = Polygon([(0.0, 0.0), (0.0, 0.04), (0.04, 0.04), (0.04, 0.0)])
    basin_gdf = gpd.GeoDataFrame({"name": ["synthetic_basin"]}, geometry=[polygon], crs="EPSG:4326")
    basin_gdf.to_file(inputs_root / "basin.shp")

    transform = from_origin(0.0, 0.04, 0.01, 0.01)
    dem = np.array(
        [
            [4300, 4400, 4500, 4600],
            [4400, 4500, 4600, 4700],
            [4500, 4600, 4700, 4800],
            [4600, 4700, 4800, 4900],
        ],
        dtype=np.float32,
    )
    write_raster(inputs_root / "dem.tif", dem, transform, nodata=-9999.0)
    write_raster(alt_inputs_root / "ALT_2006.tif", np.full((4, 4), 2.0, dtype=np.float32), transform, nodata=-9999.0)

    dates = pd.date_range("2006-01-01", "2006-12-31", freq="D")
    discharge = pd.DataFrame(
        {
            "date": dates,
            "discharge": np.linspace(5.0, 15.0, len(dates)),
        }
    )
    discharge.to_csv(inputs_root / "discharge.csv", index=False)

    gis_dir = workspace_root / "data" / "gis"
    aligned_root = workspace_root / "data" / "aligned_masked"
    alt_basin_dir = workspace_root / "data" / "alt" / "basin"

    gis_dir.mkdir(parents=True, exist_ok=True)
    (workspace_root / "results" / "runs").mkdir(parents=True, exist_ok=True)
    (workspace_root / "results" / "original_hbv").mkdir(parents=True, exist_ok=True)
    (workspace_root / "results" / "logs").mkdir(parents=True, exist_ok=True)
    (workspace_root / "results" / "cache").mkdir(parents=True, exist_ok=True)

    write_raster(gis_dir / "dem_1km.tif", dem, transform, nodata=-9999.0)
    write_raster(gis_dir / "flow_direction.tif", np.full((4, 4), 1.0, dtype=np.float32), transform, nodata=-9999.0)
    write_raster(gis_dir / "flow_accumulation.tif", np.arange(1, 17, dtype=np.float32).reshape(4, 4), transform, nodata=-9999.0)
    write_raster(gis_dir / "flow_direction_masked.tif", np.full((4, 4), 1.0, dtype=np.float32), transform, nodata=-9999.0)
    write_raster(gis_dir / "flow_accumulation_masked.tif", np.arange(1, 17, dtype=np.float32).reshape(4, 4), transform, nodata=-9999.0)
    write_uint_raster(gis_dir / "elevation_zone_low.tif", (dem < 4500).astype(np.uint8), transform)
    write_uint_raster(gis_dir / "elevation_zone_mid.tif", ((dem >= 4500) & (dem < 5000)).astype(np.uint8), transform)
    write_uint_raster(gis_dir / "elevation_zone_high.tif", (dem >= 5000).astype(np.uint8), transform)
    write_uint_raster(gis_dir / "glacier_mask.tif", np.zeros((4, 4), dtype=np.uint8), transform)

    alt_basin_dir.mkdir(parents=True, exist_ok=True)
    write_raster(alt_basin_dir / "ALT_2006.tif", np.full((4, 4), 2.0, dtype=np.float32), transform, nodata=-9999.0)

    for subdir in ["prec", "temp", "evap"]:
        (aligned_root / subdir).mkdir(parents=True, exist_ok=True)

    for index, date in enumerate(dates):
        date_str = date.strftime("%Y.%m.%d")
        prec = np.full((4, 4), 3.0 + index * 0.1, dtype=np.float32)
        temp = np.full((4, 4), -2.0 + index * 0.2, dtype=np.float32)
        evap = np.full((4, 4), 1.0, dtype=np.float32)
        write_raster(aligned_root / "prec" / f"{index:03d}_PREC_{date_str}.tif", prec, transform, nodata=-9999.0)
        write_raster(aligned_root / "temp" / f"{index:03d}_TEMP_{date_str}.tif", temp, transform, nodata=-9999.0)
        write_raster(aligned_root / "evap" / f"{index:03d}_EVAP_{date_str}.tif", evap, transform, nodata=-9999.0)

    print(f"Synthetic example generated under: {repo_root / 'examples' / 'synthetic_basin'}")


if __name__ == "__main__":
    main()
