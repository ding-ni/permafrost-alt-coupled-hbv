# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import (
    basin_paths,
    build_workspace_paths,
    default_config_path,
    ensure_workspace_dirs,
    mask_raster_to_basin,
    read_config,
)


def main():
    parser = argparse.ArgumentParser(description="Mask flow direction and flow accumulation rasters to the basin.")
    parser.add_argument("--config", default=str(default_config_path()))
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    basin = basin_paths(config)
    ensure_workspace_dirs(paths)

    mask_raster_to_basin(
        paths["gis_dir"] / "flow_direction.tif",
        paths["gis_dir"] / "flow_direction_masked.tif",
        basin["basin_shp"],
    )
    mask_raster_to_basin(
        paths["gis_dir"] / "flow_accumulation.tif",
        paths["gis_dir"] / "flow_accumulation_masked.tif",
        basin["basin_shp"],
    )
    print("Masked flow_direction.tif and flow_accumulation.tif.")


if __name__ == "__main__":
    main()
