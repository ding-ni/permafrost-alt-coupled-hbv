# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import (
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
    parser = argparse.ArgumentParser(description="Build flow direction and flow accumulation rasters.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--method", choices=["whitebox", "arcpy"], default="whitebox")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    ensure_workspace_dirs(paths)
    print_config_summary(config, paths)

    module = load_legacy_module(legacy_file("preprocess", "generate_flow_direction.py"))
    patch_module(
        module,
        {
            "GIS_DIR": str(paths["gis_dir"]),
            "DEM_INPUT": str(paths["gis_dir"] / "dem_1km.tif"),
            "DEM_FILLED": str(paths["gis_dir"] / "dem_filled.tif"),
            "FLOW_DIR": str(paths["gis_dir"] / "flow_direction.tif"),
            "FLOW_ACC": str(paths["gis_dir"] / "flow_accumulation.tif"),
        },
    )

    if args.method == "whitebox":
        module.run_whitebox()
    else:
        module.run_arcpy()
    module.verify_results()


if __name__ == "__main__":
    main()
