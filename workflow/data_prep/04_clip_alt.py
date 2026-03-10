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
    year_range,
)


def main():
    parser = argparse.ArgumentParser(description="Clip annual ALT rasters into the basin workspace.")
    parser.add_argument("--config", default=str(default_config_path()))
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    basin = basin_paths(config)
    ensure_workspace_dirs(paths)

    module = load_legacy_module(legacy_file("preprocess", "clip_alt_data.py"))
    patch_module(
        module,
        {
            "PROJECT_ROOT": str(paths["workspace_root"]),
            "BASIN_SHP": str(basin["basin_shp"]),
            "ALT_INPUT_DIR": str(basin["alt_input_dir"]),
            "ALT_OUTPUT_DIR": str(paths["alt_basin_dir"]),
            "YEARS": year_range(config),
        },
    )
    module.clip_alt_to_basin()


if __name__ == "__main__":
    main()
