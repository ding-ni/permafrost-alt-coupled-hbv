# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import (
    bbox_dict,
    build_workspace_paths,
    default_config_path,
    ensure_workspace_dirs,
    legacy_file,
    load_legacy_module,
    patch_module,
    read_config,
    temporary_argv,
    year_range,
)


def main():
    parser = argparse.ArgumentParser(description="Process precipitation rasters from MSWEP or CMFD.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--prec-source", choices=["mswep", "cmfd"], default="mswep")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    ensure_workspace_dirs(paths)
    years = list(year_range(config))

    if args.prec_source == "mswep":
        module = load_legacy_module(legacy_file("preprocess", "process_mswep.py"))
        patch_module(
            module,
            {
                "PROJECT_ROOT": str(paths["workspace_root"]),
                "MSWEP_INPUT_DIR": str(paths["raw_prec_mswep_dir"]),
                "PREC_OUTPUT_DIR": str(paths["raw_prec_daily_dir"]),
                "BASIN_BBOX": bbox_dict(config),
                "TUOTUOHE_BBOX": bbox_dict(config),
                "START_YEAR": years[0],
                "END_YEAR": years[-1],
            },
        )
        module.main()
        return

    flow_acc_masked = paths["gis_dir"] / "flow_accumulation_masked.tif"
    if not flow_acc_masked.exists():
        raise FileNotFoundError(
            "CMFD preprocessing requires flow_accumulation_masked.tif. "
            "Run workflow/data_prep/12_mask_flow_rasters.py first."
        )

    module = load_legacy_module(legacy_file("preprocess", "process_cmfd.py"))
    patch_module(
        module,
        {
            "PROJECT_ROOT": str(paths["workspace_root"]),
            "CMFD_DIR": str(paths["raw_prec_cmfd_dir"]),
            "OUT_DIR": str(paths["raw_prec_cmfd_daily_dir"]),
            "FLOW_ACC_PATH": str(flow_acc_masked),
        },
    )

    argv = ["07_process_precipitation.py"]
    if args.overwrite:
        argv.append("--overwrite")
    with temporary_argv(argv):
        module.main()


if __name__ == "__main__":
    main()
