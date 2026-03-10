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
    read_config,
    temporary_argv,
)


def main():
    parser = argparse.ArgumentParser(description="Generate glacier melt rasters from temperature and radiation.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    ensure_workspace_dirs(paths)

    module = load_legacy_module(legacy_file("glacier", "generate_glacier_melt_from_era5.py"))
    patch_module(
        module,
        {
            "PROJECT_ROOT": str(paths["workspace_root"]),
            "TEMP_DIR": str(paths["aligned_temp_dir"]),
            "RAD_DIR": str(paths["raw_solar_dir"]),
            "OUT_DIR": str(paths["glacier_melt_dir"]),
            "DEM_PATH": str(paths["gis_dir"] / "dem_1km.tif"),
            "FLOW_ACC_PATH": str(paths["gis_dir"] / "flow_accumulation_masked.tif"),
            "GLACIER_MASK_PATH": str(paths["gis_dir"] / "glacier_mask.tif"),
        },
    )

    argv = ["10_build_glacier_melt.py"]
    argv.append("--overwrite" if args.overwrite else "--no-overwrite")
    with temporary_argv(argv):
        module.main()


if __name__ == "__main__":
    main()
