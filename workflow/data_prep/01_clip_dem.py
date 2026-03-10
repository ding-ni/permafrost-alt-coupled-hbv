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
    print_config_summary,
    read_config,
)


def main():
    parser = argparse.ArgumentParser(description="Clip a DEM into the current basin workspace.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--resample-to-1km", action="store_true")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    basin = basin_paths(config)
    ensure_workspace_dirs(paths)
    print_config_summary(config, paths)

    module = load_legacy_module(legacy_file("preprocess", "clip_dem.py"))
    patch_module(
        module,
        {
            "PROJECT_ROOT": str(paths["workspace_root"]),
            "TP_DEM_FILE": str(basin["raw_dem"]),
            "BASIN_SHP": str(basin["basin_shp"]),
            "OUTPUT_DEM": str(paths["gis_dir"] / "dem_1km.tif"),
        },
    )

    ok = module.clip_dem_to_basin()
    if ok and args.resample_to_1km and hasattr(module, "resample_to_1km"):
        module.resample_to_1km()


if __name__ == "__main__":
    main()
