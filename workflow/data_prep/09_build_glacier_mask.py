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
    legacy_file,
    load_legacy_module,
    patch_module,
    read_config,
)


def main():
    parser = argparse.ArgumentParser(description="Build glacier_mask.tif from a glacier shapefile.")
    parser.add_argument("--config", default=str(default_config_path()))
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    basin = basin_paths(config)
    ensure_workspace_dirs(paths)

    if not basin["glacier_shp"]:
        raise ValueError("glacier_shp is empty in the config file.")

    module = load_legacy_module(legacy_file("glacier", "make_glacier_mask_from_shp.py"))
    patch_module(
        module,
        {
            "PROJECT_ROOT": str(paths["workspace_root"]),
            "DEM_FILE": str(paths["gis_dir"] / "dem_1km.tif"),
            "BASIN_SHP": str(basin["basin_shp"]),
            "GLACIER_SHP": str(basin["glacier_shp"]),
            "OUT_FILE": str(paths["gis_dir"] / "glacier_mask.tif"),
        },
    )
    module.main()


if __name__ == "__main__":
    main()
