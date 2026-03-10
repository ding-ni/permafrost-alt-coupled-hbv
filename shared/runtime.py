# -*- coding: utf-8 -*-
import contextlib
import importlib.util
import json
import os
import re
import sys
import time
from pathlib import Path

import geopandas as gpd
import rasterio
from rasterio.mask import mask


REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = REPO_ROOT / "config"
DEFAULT_CONFIG_PATH = CONFIG_DIR / "basin.template.json"
DEFAULT_OBJECTIVE_MONTHS = [5, 6, 7, 8, 9, 10]


def default_config_path():
    return DEFAULT_CONFIG_PATH


def resolve_path(value, base=None):
    if value is None or value == "":
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    anchor = Path(base) if base else REPO_ROOT
    return Path(os.path.abspath(str(anchor / path)))


def read_config(config_path):
    path = resolve_path(config_path)
    with path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    config["_config_path"] = str(path)
    return config


def config_base_dir(config):
    return Path(config["_config_path"]).parent


def workspace_root(config):
    return resolve_path(config["workspace_root"], base=config_base_dir(config))


def build_workspace_paths(config):
    root = workspace_root(config)
    data_root = root / "data"
    raw_root = data_root / "raw"
    raw_prec_root = raw_root / "precipitation"
    return {
        "workspace_root": root,
        "data_root": data_root,
        "gis_dir": data_root / "gis",
        "raw_root": raw_root,
        "raw_prec_root": raw_prec_root,
        "raw_prec_mswep_dir": raw_prec_root / "mswep",
        "raw_prec_daily_dir": raw_prec_root / "daily",
        "raw_prec_cmfd_dir": raw_prec_root / "CMFD",
        "raw_prec_cmfd_daily_dir": raw_prec_root / "cmfd_daily",
        "raw_temp_dir": raw_root / "temperature",
        "raw_temp_daily_dir": raw_root / "temperature" / "daily",
        "raw_evap_dir": raw_root / "evaporation",
        "raw_evap_daily_dir": raw_root / "evaporation" / "daily",
        "raw_solar_dir": raw_root / "solar_radiation",
        "raw_wind_dir": raw_root / "wind",
        "raw_dewpoint_dir": raw_root / "dewpoint",
        "aligned_dir": data_root / "aligned_masked",
        "aligned_prec_dir": data_root / "aligned_masked" / "prec",
        "aligned_prec_cmfd_dir": data_root / "aligned_masked" / "prec_cmfd",
        "aligned_temp_dir": data_root / "aligned_masked" / "temp",
        "aligned_evap_dir": data_root / "aligned_masked" / "evap",
        "glacier_melt_dir": data_root / "aligned_masked" / "glacier_melt",
        "alt_root": data_root / "alt",
        "alt_basin_dir": data_root / "alt" / "basin",
        "alt_scale_dir": data_root / "alt" / "scale_factors",
        "parameter_root": data_root / "parameters",
        "global_param_dir": data_root / "parameters" / "global",
        "alt_param_dir": data_root / "parameters" / "alt_adjusted_v6",
        "observed_dir": data_root / "observed",
        "results_root": root / "results",
        "original_results_dir": root / "results" / "original_hbv",
        "runs_dir": root / "results" / "runs",
        "figures_dir": root / "results" / "figures",
        "logs_dir": root / "results" / "logs",
        "cache_dir": root / "results" / "cache",
    }


def ensure_workspace_dirs(paths):
    for value in paths.values():
        if isinstance(value, Path):
            value.mkdir(parents=True, exist_ok=True)


def time_settings(config):
    return config["time"]


def year_range(config):
    settings = time_settings(config)
    return range(int(settings["start_year"]), int(settings["end_year"]) + 1)


def bbox_list(config):
    bbox = config["bbox"]
    return [bbox["north"], bbox["west"], bbox["south"], bbox["east"]]


def bbox_dict(config):
    bbox = config["bbox"]
    return {
        "lon_min": bbox["west"],
        "lon_max": bbox["east"],
        "lat_min": bbox["south"],
        "lat_max": bbox["north"],
    }


def basin_paths(config):
    base = config_base_dir(config)
    return {
        "basin_shp": resolve_path(config["basin_shp"], base=base),
        "raw_dem": resolve_path(config["dem_tif"], base=base),
        "obs_csv": resolve_path(config["observed_discharge_csv"], base=base),
        "glacier_shp": resolve_path(config.get("glacier_shp", ""), base=base),
        "alt_input_dir": resolve_path(config["alt_input_dir"], base=base),
    }


def expected_days(config):
    import pandas as pd

    settings = time_settings(config)
    start_date = f"{int(settings['start_year'])}-01-01"
    end_date = f"{int(settings['end_year'])}-12-31"
    return len(pd.date_range(start_date, end_date, freq="D"))


def mask_raster_to_basin(input_file, output_file, basin_shp):
    basin = gpd.read_file(basin_shp)
    with rasterio.open(input_file) as src:
        basin_reproj = basin.to_crs(src.crs)
        nodata = src.nodata if src.nodata is not None else -9999
        out_image, _ = mask(src, basin_reproj.geometry, crop=False, nodata=nodata)
        meta = src.meta.copy()
        meta.update(nodata=nodata, compress="lzw")
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with rasterio.open(output_file, "w", **meta) as dst:
            dst.write(out_image)


def legacy_file(*parts):
    return REPO_ROOT.joinpath("legacy_core", *parts)


def load_legacy_module(script_path):
    script_path = Path(script_path)
    module_name = f"legacy_{script_path.stem}_{int(time.time() * 1000)}"
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load legacy script: {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def patch_module(module, mapping):
    for key, value in mapping.items():
        setattr(module, key, value)


@contextlib.contextmanager
def temporary_argv(argv):
    old_argv = sys.argv[:]
    sys.argv = argv[:]
    try:
        yield
    finally:
        sys.argv = old_argv


def parse_params_from_result(result_path):
    result_path = resolve_path(result_path)
    if result_path.is_dir():
        metadata = result_path / "metadata.json"
        params_txt = result_path / "parameters.txt"
        if metadata.exists():
            return parse_params_from_result(metadata)
        if params_txt.exists():
            return parse_params_from_result(params_txt)
        raise FileNotFoundError(f"Run directory has no metadata.json or parameters.txt: {result_path}")

    if result_path.suffix.lower() == ".json":
        with result_path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        params = data.get("optimized_params") or data.get("params") or {}
        return {str(key).upper(): float(value) for key, value in params.items()}

    text = result_path.read_text(encoding="utf-8", errors="ignore")
    matches = re.findall(r"^([A-Za-z0-9_]+)\s*=\s*([+-]?\d+(?:\.\d+)?)$", text, flags=re.MULTILINE)
    return {key.upper(): float(value) for key, value in matches}


def load_partial_v6_params(config, result_path=None):
    if result_path:
        params = parse_params_from_result(result_path)
    else:
        params = {
            str(key).upper(): float(value)
            for key, value in config.get("improved_hbv_baseline_params", {}).items()
        }

    required = ["FC", "PERC", "K2"]
    missing = [key for key in required if key not in params]
    if missing:
        raise KeyError(f"Missing required improved HBV parameters: {missing}")
    if "K" not in params and "K0" in params:
        params["K"] = params["K0"]
    return params


def load_full_v6_params(config, result_path=None):
    source = result_path or config.get("improved_hbv_run_path", "")
    if not source:
        raise ValueError("ALT-coupled HBV requires a metadata.json file or run directory from improved HBV.")

    params = parse_params_from_result(source)
    required = [
        "TT",
        "FC",
        "BETA",
        "LP",
        "RFCF",
        "SFCF",
        "CFR",
        "CWH",
        "CFMAX_MID",
        "CFMAX_HIGH",
        "K",
        "K1",
        "K2",
        "PERC",
        "K_MUSK",
        "ICE_FACTOR",
    ]
    missing = [key for key in required if key not in params]
    if missing:
        raise KeyError(f"Improved HBV output is missing required parameters: {missing}")
    params["K0"] = params.get("K0", params["K"])
    return params


def observation_mode_for_legacy(config):
    mode = str(config.get("observation_mode", "full_year")).strip().lower()
    if mode == "full_year":
        return "full_year"
    if mode != "seasonal_window":
        raise ValueError("observation_mode must be 'full_year' or 'seasonal_window'.")

    months = sorted(int(month) for month in config.get("objective_window_months", DEFAULT_OBJECTIVE_MONTHS))
    if months != DEFAULT_OBJECTIVE_MONTHS:
        raise NotImplementedError(
            "The current public engine supports seasonal_window only for May-October "
            "(objective_window_months = [5, 6, 7, 8, 9, 10])."
        )
    return "may_oct"


def coupling_start_year(config):
    return int(config.get("coupling_start_year") or time_settings(config)["validation_start"][:4])


def print_config_summary(config, paths):
    print("=" * 72)
    print("Public Basin Configuration")
    print("=" * 72)
    print(f"Basin (EN): {config['basin_name_en']}")
    print(f"Basin (ZH): {config.get('basin_name_zh', '')}")
    print(f"Workspace: {paths['workspace_root']}")
    print(f"Catchment area (km2): {config['catchment_area_km2']}")
    print(f"Default precipitation source: {config.get('default_precip_source', 'mswep')}")
    print(f"Observation mode: {config.get('observation_mode', 'full_year')}")
    print("=" * 72)
