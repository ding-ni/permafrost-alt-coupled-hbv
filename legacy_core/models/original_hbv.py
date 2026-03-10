# -*- coding: utf-8 -*-
"""
原始分布式HBV模型对比实验 V2

与改进型分布式HBV的区别：
1. 禁用冰川模块（无ICE_FACTOR）
2. 使用单一CFMAX（不分区mid/high）
3. UZL和X_MUSK作为率定参数（改进型固定）
4. e_corr固定为0（与改进型相同，因为都用FAO-56 PET）
5. 目标函数：简单的率定期NSE（改进型用多目标复合）
6. 无径流来源追踪

运行方式:
  python 00_run_original_hbv_v2.py --maxiter 80 --popsize 12

输出:
  结果/原始HBV对比/original_hbv_summary.json
  结果/原始HBV对比/original_hbv_timeseries.csv
"""

import sys
import os
import time
import json
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from glob import glob
from math import sin, radians
import rasterio
from scipy.optimize import differential_evolution
from numba import njit
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 配置
# =============================================================================
PROJECT_ROOT = Path(".").resolve()
THESIS_ROOT = PROJECT_ROOT

PREC_DIR = PROJECT_ROOT / 'data' / 'aligned_masked' / 'prec'
TEMP_DIR = PROJECT_ROOT / 'data' / 'aligned_masked' / 'temp'
EVAP_DIR = PROJECT_ROOT / 'data' / 'aligned_masked' / 'evap'
OBS_FILE = PROJECT_ROOT / 'data' / 'observed' / 'discharge_example.csv'
FLOW_ACC_PATH = PROJECT_ROOT / 'data' / 'gis' / 'flow_accumulation_masked.tif'

RESULTS_DIR = THESIS_ROOT / '结果' / '原始HBV对比'

CATCHMENT_AREA = 15924  # km2
EPS = 1e-12
MUSK_DT = 1.0

# 时间配置
WARMUP_START = '2006-01-01'
WARMUP_END = '2008-12-31'
CALIB_START = '2009-01-01'
CALIB_END = '2017-12-31'
VALID_START = '2018-01-01'
VALID_END = '2020-12-31'
SIM_START = CALIB_START
SIM_END = VALID_END
OBS_MODE_OVERRIDE = 'auto'  # auto | full_year | may_oct

# 目标函数：简单NSE（原始HBV标准做法）
# 不使用多目标复合、不使用验证期数据、不使用季节性惩罚

# =============================================================================
# 原始HBV参数边界（16参数，e_corr固定为0）
# 与改进型HBV的区别：单一CFMAX、无冰川模块、UZL和X_MUSK率定
# e_corr固定为0（因为使用FAO-56 PET，无需温度修正）
# =============================================================================
PARAM_NAMES = [
    'tt', 'rfcf', 'sfcf', 'cfmax', 'cwh', 'cfr',
    'fc', 'beta', 'lp',
    'k', 'k1', 'k2', 'uzl', 'perc', 'k_musk', 'x_musk'
]

PARAM_BOUNDS = [
    (-3.0, 3.0),      # tt
    (0.8, 1.5),       # rfcf
    (0.5, 2.0),       # sfcf
    (1.0, 10.0),      # cfmax
    (0.01, 0.2),      # cwh
    (0.001, 0.1),     # cfr
    (50.0, 800.0),    # fc
    (0.5, 5.0),       # beta
    (0.2, 1.0),       # lp
    (0.02, 0.5),      # k
    (0.01, 0.2),      # k1
    (0.0005, 0.02),   # k2
    (10.0, 100.0),    # uzl
    (0.1, 5.0),       # perc
    (0.625, 2.5),     # k_musk
    (0.1, 0.4),       # x_musk
]

# 固定参数
E_CORR = 0.0  # FAO-56 PET无需温度修正

INIT_ST = np.array([0.0, 5.0, 0.0, 0.0, 0.0], dtype=np.float64)

# =============================================================================
# 全局数据
# =============================================================================
PREC_3D = None
TEMP_3D = None
ET_3D = None
LL_TEMP_3D = None
FLOW_ACC = None
VALID_CELLS = None
CELL_SCALE = None
Q_OBS = None
Q_OBS_OBJ = None
WARMUP_DAYS = None
SIM_DATES = None
CALIB_MASK = None
VALID_MASK = None
OBS_EVAL_MASK = None
OBS_EVAL_DESC = None

# =============================================================================
# Numba加速函数
# =============================================================================
@njit(cache=True, fastmath=True, nogil=True)
def hbv_cell_original(prec, temp, et, ll_temp, par, init_st, e_corr):
    """
    原始HBV单栅格模拟（无冰川模块，无径流追踪）
    16个率定参数 + e_corr固定为0
    """
    n = len(prec)
    q_uz = np.zeros(n, dtype=np.float32)
    q_lz = np.zeros(n, dtype=np.float32)

    tt, rfcf, sfcf, cfmax, cwh, cfr = par[0], par[1], par[2], par[3], par[4], par[5]
    fc, beta, lp = par[6], par[7], par[8]
    k, k1, k2, uzl, perc = par[9], par[10], par[11], par[12], par[13]

    fc = max(fc, 10.0)
    beta = max(min(beta, 10.0), 0.1)
    lp = max(min(lp, 0.99), 0.01)

    sp, sm, uz, lz, wc = init_st[0], init_st[1], init_st[2], init_st[3], init_st[4]

    for i in range(n):
        p = prec[i]
        t = temp[i]
        e = et[i]
        tm = ll_temp[i]

        if np.isnan(p) or np.isnan(t) or np.isnan(e):
            q_uz[i] = 0.0
            q_lz[i] = 0.0
            continue

        p = max(p, 0.0)
        e = max(e, 0.0)

        # 降水相态
        if t <= tt:
            rf, sf = 0.0, p * sfcf
        else:
            rf, sf = p * rfcf, 0.0

        # 积雪模块
        if t > tt:
            avail_snow = sp + sf
            melt = min(cfmax * (t - tt), avail_snow)
            sp = max(avail_snow - melt, 0.0)
            wc_int = wc + melt + rf
        else:
            refr = min(cfr * cfmax * (tt - t), wc + rf)
            sp = sp + sf + refr
            wc_int = max(wc - refr + rf, 0.0)

        sp = min(sp, 10000.0)

        if wc_int > cwh * sp:
            inf = wc_int - cwh * sp
            wc = cwh * sp
        else:
            inf = 0.0
            wc = wc_int

        # 土壤模块（带e_corr温度修正）
        sm_ratio = min(max(sm / fc, 0.0), 1.0)
        r = (sm_ratio ** beta) * inf

        ep_adj = max((1.0 + (t - tm) * e_corr) * e, 0.0)
        lp_fc = lp * fc
        if lp_fc > 0.001:
            ea = min(ep_adj, (sm / lp_fc) * ep_adj)
        else:
            ea = ep_adj
        ea = min(ea, sm)

        uz_int = uz + r
        sm = max(min(sm + inf - r - ea, fc), 0.0)

        # 响应函数
        perc_actual = min(perc, uz_int)
        uz_int2 = max(uz_int - perc_actual, 0.0)

        q0 = k * max(uz_int2 - uzl, 0.0)
        q1 = k1 * uz_int2

        if q0 + q1 > uz_int2:
            q0 = uz_int2 * 0.67
            q1 = uz_int2 * 0.33

        uz = max(uz_int2 - (q0 + q1), 0.0)

        lz_int = lz + perc_actual
        q2 = k2 * lz_int
        if q2 > lz_int:
            q2 = lz_int
        lz = max(lz_int - q2, 0.0)

        uz = min(uz, 10000.0)
        lz = min(lz, 10000.0)

        q_uz[i] = q0 + q1
        q_lz[i] = q2

    return q_uz, q_lz


@njit(cache=True, fastmath=True, nogil=True)
def run_all_cells_original(prec_3d, temp_3d, et_3d, ll_temp_3d, par, init_st,
                           valid_cells, cell_scale, e_corr):
    """分布式运行所有栅格"""
    n_cells = valid_cells.shape[0]
    ts = prec_3d.shape[2]
    q_total = np.zeros(ts, dtype=np.float64)

    for idx in range(n_cells):
        x = valid_cells[idx, 0]
        y = valid_cells[idx, 1]
        scale = cell_scale[idx]

        prec = prec_3d[x, y, :]
        temp = temp_3d[x, y, :]
        et = et_3d[x, y, :]
        ll_temp = ll_temp_3d[x, y, :]

        q_uz, q_lz = hbv_cell_original(prec, temp, et, ll_temp, par, init_st, e_corr)

        for t in range(ts):
            q_total[t] += (q_uz[t] + q_lz[t]) * scale

    return q_total


@njit(cache=True, fastmath=True, nogil=True)
def muskingum_route(q_in, k_musk, x_musk, dt):
    """Muskingum汇流"""
    n = len(q_in)
    q_out = np.zeros(n, dtype=np.float64)

    denom = 2.0 * k_musk * (1.0 - x_musk) + dt
    c0 = (dt - 2.0 * k_musk * x_musk) / denom
    c1 = (dt + 2.0 * k_musk * x_musk) / denom
    c2 = (2.0 * k_musk * (1.0 - x_musk) - dt) / denom

    c0 = max(c0, 0.0)
    c1 = max(c1, 0.0)
    c2 = max(c2, 0.0)

    c_sum = c0 + c1 + c2
    if c_sum > 0:
        c0 /= c_sum
        c1 /= c_sum
        c2 /= c_sum

    q_out[0] = q_in[0]
    for i in range(1, n):
        q_out[i] = c0 * q_in[i] + c1 * q_in[i-1] + c2 * q_out[i-1]

    return q_out


# =============================================================================
# 数据加载
# =============================================================================
def compute_row_areas_km2(transform, height):
    pixel_width_deg = transform.a
    pixel_height_deg = -transform.e
    top_lat = transform.f
    r = 6371000.0
    dlon = radians(pixel_width_deg)
    row_areas = np.zeros(height, dtype=np.float64)
    for row in range(height):
        lat_n = top_lat - row * pixel_height_deg
        lat_s = top_lat - (row + 1) * pixel_height_deg
        row_areas[row] = (r ** 2) * dlon * (sin(radians(lat_n)) - sin(radians(lat_s))) / 1e6
    return row_areas


def load_raster_stack(directory, start_date, end_date):
    dates = pd.date_range(start_date, end_date, freq='D')
    all_files = glob(str(directory / '*.tif'))
    file_date_map = {}
    for f in all_files:
        fname = os.path.basename(f)
        parts = fname.replace('.tif', '').split('_')
        date_str = parts[-1] if len(parts) >= 3 else parts[0]
        try:
            fdate = pd.to_datetime(date_str, format='%Y.%m.%d')
            file_date_map[fdate] = f
        except:
            continue

    first_file = next(iter(file_date_map.values()))
    with rasterio.open(first_file) as src:
        rows, cols = src.height, src.width
        transform = src.transform

    data = np.full((rows, cols, len(dates)), np.nan, dtype=np.float32)
    for i, d in enumerate(dates):
        if d in file_date_map:
            with rasterio.open(file_date_map[d]) as src:
                arr = src.read(1).astype(np.float32)
                nodata = src.nodata
                if nodata is not None:
                    arr[arr == nodata] = np.nan
                arr[arr < -9000] = np.nan
                arr[arr > 1e10] = np.nan
                data[:, :, i] = arr
    return data, transform


def detect_obs_eval_mask(q_obs, sim_dates, mode_override='auto'):
    full_mask = np.ones(len(sim_dates), dtype=bool)
    may_oct_mask = (sim_dates.month >= 5) & (sim_dates.month <= 10)
    mode = str(mode_override or 'auto').strip().lower()

    if mode == 'full_year':
        return full_mask, 'full_year_forced'
    if mode == 'may_oct':
        return may_oct_mask, 'mayoct_forced'
    if mode != 'auto':
        print(f"[WARN] Unknown OBS_MODE_OVERRIDE={mode_override!r}, fallback to auto.")

    try:
        df_frac = pd.DataFrame({'date': sim_dates, 'q': q_obs})
        df_frac['year'] = df_frac['date'].dt.year
        df_frac['month'] = df_frac['date'].dt.month
        frac_list = []
        for _, g in df_frac.groupby('year'):
            qy = g['q'].values.astype(np.float64)
            tot = np.nansum(qy)
            if tot <= 0:
                continue
            mayoct = np.nansum(g.loc[(g['month'] >= 5) & (g['month'] <= 10), 'q'].values.astype(np.float64))
            frac_list.append(float(mayoct / tot))
        frac_arr = np.asarray(frac_list, dtype=np.float64)
        frac_mean = float(np.nanmean(frac_arr)) if frac_arr.size else float("nan")
        frac_std = float(np.nanstd(frac_arr)) if frac_arr.size else float("nan")
        if np.isfinite(frac_mean) and np.isfinite(frac_std) and (frac_std < 1e-4) and (abs(frac_mean - 0.967) < 0.01):
            return may_oct_mask, f"mayoct_only_yearbook(frac={frac_mean:.3f}±{frac_std:.3g})"
    except Exception as exc:
        print(f"[WARN] Obs mode detection failed, fallback to full-year eval: {exc}")

    return full_mask, 'full_year'


def load_all_data():
    global PREC_3D, TEMP_3D, ET_3D, LL_TEMP_3D, FLOW_ACC
    global VALID_CELLS, CELL_SCALE, Q_OBS, Q_OBS_OBJ, WARMUP_DAYS, SIM_DATES
    global CALIB_MASK, VALID_MASK, OBS_EVAL_MASK, OBS_EVAL_DESC

    print("=" * 70)
    print("Loading data...")
    print("=" * 70)

    print("[1/5] Loading basin mask...")
    with rasterio.open(FLOW_ACC_PATH) as src:
        flow_acc = src.read(1).astype(np.float32)
        nodata = src.nodata
        if nodata is not None:
            flow_acc[flow_acc == nodata] = np.nan
        transform = src.transform
        rows, cols = src.height, src.width

    valid_mask = ~np.isnan(flow_acc)
    FLOW_ACC = flow_acc

    valid_idx = np.where(valid_mask)
    VALID_CELLS = np.column_stack((valid_idx[0], valid_idx[1])).astype(np.int32)
    n_valid = len(VALID_CELLS)
    print(f"      Valid cells: {n_valid}")

    row_areas = compute_row_areas_km2(transform, rows)
    px_area = np.zeros((rows, cols), dtype=np.float64)
    for i in range(rows):
        px_area[i, :] = row_areas[i]

    px_tot_area = np.sum(px_area[valid_mask])
    area_coef = CATCHMENT_AREA / px_tot_area
    conversion_factor = 86.4
    CELL_SCALE = np.zeros(n_valid, dtype=np.float64)
    for idx in range(n_valid):
        x, y = VALID_CELLS[idx]
        CELL_SCALE[idx] = px_area[x, y] * area_coef / conversion_factor

    print("[2/5] Loading precipitation...")
    PREC_3D, _ = load_raster_stack(PREC_DIR, WARMUP_START, SIM_END)
    PREC_3D = np.nan_to_num(PREC_3D, nan=0.0)
    print(f"      Mean: {np.mean(PREC_3D[valid_mask]):.2f} mm/d")

    print("[3/5] Loading temperature...")
    TEMP_3D, _ = load_raster_stack(TEMP_DIR, WARMUP_START, SIM_END)
    ll_temp = np.nanmean(TEMP_3D, axis=2, keepdims=True)
    ll_temp = np.broadcast_to(ll_temp, TEMP_3D.shape).copy()
    LL_TEMP_3D = ll_temp
    TEMP_3D = np.where(np.isnan(TEMP_3D), ll_temp, TEMP_3D)
    print(f"      Mean: {np.nanmean(TEMP_3D[valid_mask]):.2f} degC")

    print("[4/5] Loading evapotranspiration...")
    ET_3D, _ = load_raster_stack(EVAP_DIR, WARMUP_START, SIM_END)
    ET_3D = np.nan_to_num(ET_3D, nan=0.0)
    print(f"      Mean: {np.mean(ET_3D[valid_mask]):.2f} mm/d")

    print("[5/5] Loading observed discharge...")
    all_dates = pd.date_range(WARMUP_START, SIM_END, freq='D')
    SIM_DATES = pd.date_range(SIM_START, SIM_END, freq='D')
    WARMUP_DAYS = (pd.to_datetime(SIM_START) - pd.to_datetime(WARMUP_START)).days

    if OBS_FILE.exists():
        df_obs = pd.read_csv(OBS_FILE)
        df_obs['date'] = pd.to_datetime(df_obs['date'])
        df_obs.set_index('date', inplace=True)
        Q_OBS = df_obs.reindex(SIM_DATES)['discharge (m3/s)'].values.astype(np.float64)
        print(f"      Mean: {np.nanmean(Q_OBS):.2f} m3/s")
    else:
        Q_OBS = np.zeros(len(SIM_DATES))

    CALIB_MASK = (SIM_DATES >= pd.to_datetime(CALIB_START)) & (SIM_DATES <= pd.to_datetime(CALIB_END))
    VALID_MASK = (SIM_DATES >= pd.to_datetime(VALID_START)) & (SIM_DATES <= pd.to_datetime(VALID_END))
    OBS_EVAL_MASK, OBS_EVAL_DESC = detect_obs_eval_mask(Q_OBS, SIM_DATES, OBS_MODE_OVERRIDE)
    Q_OBS_OBJ = Q_OBS.copy()
    Q_OBS_OBJ[~OBS_EVAL_MASK] = np.nan

    print(f"      Obs mode: {OBS_EVAL_DESC}")
    if OBS_EVAL_DESC.startswith('mayoct'):
        print("      Objective window: May-Oct (5-10); full-period metrics will also be exported.")
    else:
        print("      Objective window: full year.")

    print("=" * 70)
    print(f"Data loaded: {PREC_3D.shape[2]} days")
    print(f"Warmup: {WARMUP_DAYS} days, Simulation: {len(SIM_DATES)} days")
    print("=" * 70)


# =============================================================================
# 模型运行
# =============================================================================
def run_hbv(params):
    par = np.array(params[:14], dtype=np.float64)  # 14个HBV参数
    k_musk = params[14]
    x_musk = params[15]

    q_total = run_all_cells_original(
        PREC_3D, TEMP_3D, ET_3D, LL_TEMP_3D,
        par, INIT_ST, VALID_CELLS, CELL_SCALE, E_CORR
    )

    q_total = q_total[WARMUP_DAYS:]
    q_routed = muskingum_route(q_total, k_musk, x_musk, MUSK_DT)

    return q_routed


def nse(obs, sim):
    mask = np.isfinite(obs) & np.isfinite(sim)
    if mask.sum() == 0:
        return -999.0
    o, s = obs[mask], sim[mask]
    ss_err = np.sum((o - s) ** 2)
    ss_tot = np.sum((o - o.mean()) ** 2)
    return 1.0 - ss_err / ss_tot if ss_tot > 0 else -999.0


def nse_log(obs, sim):
    mask = np.isfinite(obs) & np.isfinite(sim) & (obs > 0) & (sim > 0)
    if mask.sum() == 0:
        return -999.0
    o, s = np.log1p(obs[mask]), np.log1p(sim[mask])
    ss_err = np.sum((o - s) ** 2)
    ss_tot = np.sum((o - o.mean()) ** 2)
    return 1.0 - ss_err / ss_tot if ss_tot > 0 else -999.0


def pbias(obs, sim):
    mask = np.isfinite(obs) & np.isfinite(sim)
    if mask.sum() == 0:
        return 0.0
    o, s = obs[mask], sim[mask]
    obs_sum = np.sum(o)
    return 100.0 * np.sum(s - o) / obs_sum if obs_sum > 0 else 0.0


def rmse(obs, sim):
    mask = np.isfinite(obs) & np.isfinite(sim)
    if mask.sum() == 0:
        return float("nan")
    o, s = obs[mask], sim[mask]
    return float(np.sqrt(np.mean((o - s) ** 2)))


def objective_function(params):
    """原始HBV目标函数：简单的率定期NSE"""
    try:
        q_sim = run_hbv(params)
    except Exception:
        return 1e6

    if np.any(np.isnan(q_sim)) or np.any(np.isinf(q_sim)):
        return 1e6

    # 仅使用观测目标窗口计算NSE（全年实测或 5-10 月实测窗口）
    cal_mask = CALIB_MASK & OBS_EVAL_MASK
    nse_cal = nse(Q_OBS_OBJ[cal_mask], q_sim[cal_mask])

    if nse_cal < -10:
        return 1e6

    return -nse_cal  # 最小化负NSE = 最大化NSE


# =============================================================================
# 主函数
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description='Original Distributed HBV Baseline')
    parser.add_argument('--maxiter', type=int, default=30, help='DE max iterations')
    parser.add_argument('--popsize', type=int, default=5, help='DE population size factor')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    args = parser.parse_args()

    print("\n" + "=" * 70)
    print("Original Distributed HBV Model (16 params, e_corr=0)")
    print("=" * 70)

    start_time = time.time()
    load_all_data()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # JIT warmup
    print("\nJIT warmup...")
    test_params = np.array([
        0.0, 1.0, 1.0, 4.0, 0.05, 0.05,   # tt, rfcf, sfcf, cfmax, cwh, cfr
        500.0, 2.0, 0.5,                   # fc, beta, lp
        0.1, 0.05, 0.005, 30.0, 1.0,      # k, k1, k2, uzl, perc
        1.0, 0.2                           # k_musk, x_musk
    ], dtype=np.float64)
    _ = run_hbv(test_params)
    print("JIT warmup complete")

    # Optimization
    print(f"\nStarting calibration (DE)...")
    print(f"Population: {args.popsize}, Max iterations: {args.maxiter}")
    print(f"Parameters: {len(PARAM_NAMES)}")

    result = differential_evolution(
        objective_function,
        bounds=PARAM_BOUNDS,
        seed=args.seed,
        maxiter=args.maxiter,
        popsize=args.popsize,
        mutation=(0.5, 1.0),
        recombination=0.7,
        tol=1e-6,
        polish=True,
        disp=True
    )

    best_params = result.x

    # Final run
    print("\nRunning with best parameters...")
    q_sim = run_hbv(best_params)

    cal_mask_obj = CALIB_MASK & OBS_EVAL_MASK
    val_mask_obj = VALID_MASK & OBS_EVAL_MASK
    cal_mask_full = CALIB_MASK
    val_mask_full = VALID_MASK

    metrics_objective = {
        'nse_cal': float(nse(Q_OBS_OBJ[cal_mask_obj], q_sim[cal_mask_obj])),
        'nse_val': float(nse(Q_OBS_OBJ[val_mask_obj], q_sim[val_mask_obj])),
        'pbias_cal': float(pbias(Q_OBS_OBJ[cal_mask_obj], q_sim[cal_mask_obj])),
        'pbias_val': float(pbias(Q_OBS_OBJ[val_mask_obj], q_sim[val_mask_obj])),
        'rmse_cal': rmse(Q_OBS_OBJ[cal_mask_obj], q_sim[cal_mask_obj]),
        'rmse_val': rmse(Q_OBS_OBJ[val_mask_obj], q_sim[val_mask_obj]),
    }

    metrics_full_period = {
        'nse_cal': float(nse(Q_OBS[cal_mask_full], q_sim[cal_mask_full])),
        'nse_val': float(nse(Q_OBS[val_mask_full], q_sim[val_mask_full])),
        'pbias_cal': float(pbias(Q_OBS[cal_mask_full], q_sim[cal_mask_full])),
        'pbias_val': float(pbias(Q_OBS[val_mask_full], q_sim[val_mask_full])),
        'rmse_cal': rmse(Q_OBS[cal_mask_full], q_sim[cal_mask_full]),
        'rmse_val': rmse(Q_OBS[val_mask_full], q_sim[val_mask_full]),
    }

    metrics = metrics_full_period

    print("\n" + "-" * 50)
    print("Performance Metrics (objective window):")
    print(f"  Calibration: NSE={metrics_objective['nse_cal']:.4f}, PBIAS={metrics_objective['pbias_cal']:.1f}%")
    print(f"  Validation:  NSE={metrics_objective['nse_val']:.4f}, PBIAS={metrics_objective['pbias_val']:.1f}%")
    print("Performance Metrics (full period):")
    print(f"  Calibration: NSE={metrics_full_period['nse_cal']:.4f}, PBIAS={metrics_full_period['pbias_cal']:.1f}%")
    print(f"  Validation:  NSE={metrics_full_period['nse_val']:.4f}, PBIAS={metrics_full_period['pbias_val']:.1f}%")
    print("-" * 50)

    # Save timeseries
    df = pd.DataFrame({
        'date': SIM_DATES,
        'q_sim': q_sim,
        'q_obs': Q_OBS,
        'q_obs_objective': Q_OBS_OBJ,
    })
    ts_file = RESULTS_DIR / 'original_hbv_timeseries.csv'
    df.to_csv(ts_file, index=False)
    print(f"\nTimeseries saved: {ts_file}")

    # Save summary
    summary = {
        'run_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'model': 'HBV_Distributed_Original',
        'model_type': 'distributed',
        'n_params': 16,
        'description': 'Original Hapi distributed HBV (16 params), no glacier, single CFMAX, e_corr=0 (FAO-56 PET)',
        'optimization': {
            'method': 'differential_evolution',
            'objective': 'NSE (objective window only)',
            'pop_size': args.popsize,
            'max_iter': args.maxiter,
            'seed': args.seed,
        },
        'obs_mode': {
            'override': OBS_MODE_OVERRIDE,
            'detected': OBS_EVAL_DESC,
            'objective_window': 'may_oct' if OBS_EVAL_DESC.startswith('mayoct') else 'full_year',
            'comparison_metrics': 'full_period',
        },
        'fixed_params': {'e_corr': E_CORR},
        'parameters': {name: float(val) for name, val in zip(PARAM_NAMES, best_params)},
        'metrics': metrics,
        'metrics_objective_window': metrics_objective,
        'metrics_full_period': metrics_full_period,
        'time_config': {
            'warmup': f"{WARMUP_START} to {WARMUP_END}",
            'calibration': f"{CALIB_START} to {CALIB_END}",
            'validation': f"{VALID_START} to {VALID_END}",
        }
    }

    summary_file = RESULTS_DIR / 'original_hbv_summary.json'
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"Summary saved: {summary_file}")

    elapsed = time.time() - start_time
    print(f"\nCompleted in {elapsed / 60:.1f} minutes")

    return summary


if __name__ == "__main__":
    main()

