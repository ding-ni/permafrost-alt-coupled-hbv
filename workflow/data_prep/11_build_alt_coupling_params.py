# -*- coding: utf-8 -*-
import argparse
import sys
from pathlib import Path

import rasterio

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from shared.runtime import (
    build_workspace_paths,
    default_config_path,
    ensure_workspace_dirs,
    legacy_file,
    load_legacy_module,
    load_partial_v6_params,
    patch_module,
    read_config,
    year_range,
)


def format_tag(value):
    return str(value).replace("-", "m").replace(".", "p")


def build_run_tag(formula, alpha, power, include_k, ref_start, ref_end):
    default_signature = (
        formula == "inverse"
        and abs(alpha - 2.0) < 1e-12
        and abs(power - 1.0) < 1e-12
        and not include_k
        and ref_start == 2006
        and ref_end == 2017
    )
    if default_signature:
        return ""
    parts = [formula, f"a{format_tag(alpha)}", f"p{format_tag(power)}"]
    if include_k:
        parts.append("K")
    parts.append(f"ref{ref_start}_{ref_end}")
    return "_" + "_".join(parts)


def main():
    parser = argparse.ArgumentParser(description="Build ALT reference fields, scale factors, and coupling parameters.")
    parser.add_argument("--config", default=str(default_config_path()))
    parser.add_argument("--improved-run", default="")
    parser.add_argument("--formula", choices=["inverse", "forward"], default="inverse")
    parser.add_argument("--alpha", type=float, default=2.0)
    parser.add_argument("--power", type=float, default=1.0)
    parser.add_argument("--generate-k", action="store_true")
    args = parser.parse_args()

    config = read_config(args.config)
    paths = build_workspace_paths(config)
    ensure_workspace_dirs(paths)
    params = load_partial_v6_params(config, args.improved_run or None)

    module = load_legacy_module(legacy_file("preprocess", "alt_parameter_correction.py"))
    ref_cfg = config["alt_reference"]
    years = list(year_range(config))
    run_tag = build_run_tag(
        args.formula,
        float(args.alpha),
        float(args.power),
        bool(args.generate_k),
        int(ref_cfg["start_year"]),
        int(ref_cfg["end_year"]),
    )

    patch_module(
        module,
        {
            "PROJECT_ROOT": str(paths["workspace_root"]),
            "DEM_FILE": str(paths["gis_dir"] / "dem_1km.tif"),
            "ALT_DIR": str(paths["alt_basin_dir"]),
            "SF_DIR": str(paths["alt_scale_dir"]),
            "GLOBAL_PARAM_DIR": str(paths["global_param_dir"]),
            "ALT_PARAM_DIR": str(paths["alt_param_dir"]),
            "FIGURES_DIR": str(paths["figures_dir"]),
            "FC_CALIBRATED": float(params["FC"]),
            "PERC_CALIBRATED": float(params["PERC"]),
            "K2_CALIBRATED": float(params["K2"]),
            "K_CALIBRATED": float(params.get("K", params.get("K0", 0.0))),
            "SCALE_FORMULA": args.formula,
            "AMPLIFICATION_ALPHA": float(args.alpha),
            "SCALE_POWER": float(args.power),
            "GENERATE_K": bool(args.generate_k),
            "YEARS": years,
            "REF_START_YEAR": int(ref_cfg["start_year"]),
            "REF_END_YEAR": int(ref_cfg["end_year"]),
            "RUN_TAG": run_tag,
        },
    )

    with rasterio.open(paths["gis_dir"] / "dem_1km.tif") as dem:
        dem_profile = dem.profile.copy()
        dem_shape = (dem.height, dem.width)

    alt_ref = module.calculate_alt_reference(dem_profile, dem_shape)
    if alt_ref is None:
        raise RuntimeError("ALT reference field generation failed.")
    module.calculate_scale_factors(alt_ref, dem_profile, dem_shape)
    module.generate_adjusted_parameters(alt_ref, dem_profile)
    module.visualize_results()


if __name__ == "__main__":
    main()
