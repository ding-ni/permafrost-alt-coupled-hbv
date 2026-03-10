# -*- coding: utf-8 -*-
"""
目标流域 HBV 半分布式模型（分布式驱动 + 分区参数）+ 冰川融水（inline）
V6 空间分区率定脚本

核心特点：
1. CFMAX 空间分区：中海拔区(4500-5000m)与高海拔区(>5000m)独立率定。
2. 冰川消融：在冰川像元上，积雪水当量(SWE, mm w.e.)低于阈值后才允许裸冰消融；并采用
   “先融雪、后融冰”的单一度日融化潜势，避免同一日内能量被重复计入导致冰川贡献虚高。
3. 蒸散输入：ET_3D 为基于 ERA5 的 FAO-56 Penman–Monteith 计算得到的日尺度潜在蒸散(PET)栅格，
   因此关闭 HBV 传统的温度修正项（ALPHA/ECORR 固定为 0）。

率定参数（16个）：
TT, FC, BETA, LP, RFCF, SFCF, CFR, CWH, CFMAX_mid, CFMAX_high, K, K1, K2, PERC, K_MUSK, ICE_FACTOR
"""
import sys
import argparse
import os
import time
import numpy as np
import pandas as pd
import json
from datetime import datetime, timedelta
from glob import glob
from math import sin, radians
import rasterio
from scipy.optimize import differential_evolution
from numba import njit
from multiprocessing.pool import ThreadPool
import multiprocessing as mp
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 配置
# =============================================================================
PROJECT_ROOT = r'.'
PREC_DIR_MSWEP = os.path.join(PROJECT_ROOT, 'data', 'aligned_masked', 'prec')
PREC_DIR_CMFD = os.path.join(PROJECT_ROOT, 'data', 'aligned_masked', 'prec_cmfd')
PREC_DIR = PREC_DIR_MSWEP
TEMP_DIR = os.path.join(PROJECT_ROOT, 'data', 'aligned_masked', 'temp')
EVAP_DIR = os.path.join(PROJECT_ROOT, 'data', 'aligned_masked', 'evap')
OBS_FILE = os.path.join(PROJECT_ROOT, 'data', 'observed', 'discharge_example.csv')
FLOW_ACC_PATH = os.path.join(PROJECT_ROOT, 'data', 'gis', 'flow_accumulation_masked.tif')
GLACIER_MELT_DIR = os.path.join(PROJECT_ROOT, 'data', 'aligned_masked', 'glacier_melt')
GLACIER_MASK_PATH = os.path.join(PROJECT_ROOT, 'data', 'gis', 'glacier_mask.tif')
LOG_DIR = os.path.join(PROJECT_ROOT, 'results', 'logs')
CACHE_DIR = os.path.join(PROJECT_ROOT, 'results', 'cache')

CATCHMENT_AREA = 15924  # km2
GLACIER_MELT_FACTOR = 1.0
SWE_ICE_THRESHOLD_MM = 10.0  # mm w.e.; SWE<=阈值时，冰川像元允许裸冰消融
EPS = 1e-12  # 数值稳定项，避免 0/0 引发 NaN（尤其在 fastmath 下）

# 冬季低流量/缺测期处理思路（年鉴径流序列）：
# - 若识别到 11-4 月为缺测/推算期，则目标函数仅在 5-10 月实测窗口拟合；
# - 同时通过“5-10 月占比”与“推算期月体积”软约束抑制冬季非物理高流量。

# Muskingum 路由的稳定域（保证系数非负，从而分量路由仍保持线性与闭合）
# 对于 dt=1 day, X=0.2：要求 c0>=0 且 c2>=0，可得 K ∈ [dt/(2(1-X)), dt/(2X)] = [0.625, 2.5]
MUSK_DT = 1.0
X_MUSK_FIXED = 0.2
K_MUSK_MIN_SAFE = MUSK_DT / (2.0 * (1.0 - X_MUSK_FIXED)) + 1e-6
K_MUSK_MAX_SAFE = MUSK_DT / (2.0 * X_MUSK_FIXED) - 1e-6

# 无效解/异常解的惩罚值：必须足够大，避免优化器“偏爱”无效解导致提前收敛。
BAD_OBJ = 1e6

# 目标函数（推荐且可解释）：避免对源分离占比做“强行压制”，用多指标约束水文过程形态
# - NSE(Q): 约束洪峰与总体过程
# - NSE(log1p(Q)): 约束低流量/枯季过程（对0值友好）
# - PBIAS: 总体水量偏差
# - 验证期加权: 抑制过拟合，提升跨期稳健性
OBJ_W_VAL = 0.35
OBJ_W_LOG = 0.35
OBJ_W_PBIAS_CAL = 0.10
OBJ_W_PBIAS_VAL = 0.10

# 年鉴径流序列缺测/推算期识别后，默认仅用 5-10 月实测窗口进行拟合。
# 但若完全不约束 11-4 月，优化器可能将水量“搬运”到冬季，导致不合理的冬季高流量。
# 这里引入一个弱季节性约束：要求模拟的 5-10 月径流占比不低于 (年鉴占比 - slack)。
YEARBOOK_MAYOCT_FRAC_DEFAULT = 0.967
YEARBOOK_MAYOCT_FRAC_SLACK = 0.02
OBJ_W_SEASON = 100.0

# 年鉴推算期信息的“正确用法”：按月体积（或季节占比）约束，而不是把推算日值当作真实逐日过程去拟合 NSE。
OFFSEASON_MONTHS = (11, 12, 1, 2, 3, 4)
OFFSEASON_ZERO_FRAC_SLACK = 0.001  # 对目标为 0 的月份（12-2）允许少量基流：占全年体积的阈值（每月）
# 默认关闭“推算期月体积约束”：很多情况下年鉴推算期并不具备逐月可检验性，
# 强行约束会显著牺牲 5-10 月实测期拟合。若后续需要再启用，可将权重调大（如 1-10）。
OBJ_W_OFFSEASON_MONTHLY = 0.0

# 冬季冰源约束（仅年鉴模式启用；soft penalty，避免“冬季裸冰融/冰源出流”非物理解）
WINTER_ICE_MAX_MEAN_M3S = 0.5
OBJ_W_WINTER_ICE = 0.5

# 冰川融水（GM_*.tif）独立约束（可解释、可复现）
# - 用于避免“通过调参把缺水部分强行补给到冰川融水”，导致冰川贡献/消融量超出物理合理范围。
# - 该约束不直接改变模型结构，只在目标函数中对“模拟冰源出流量级”做软惩罚。
# - 默认启用：若没有 GM 数据目录/缓存则自动降级为不启用（penalty=0）。
OBJ_W_GM_VOL = 1.0
GM_VOL_RATIO_MIN = 0.7  # 允许模拟冰源体积 / GM 体积 的下限
GM_VOL_RATIO_MAX = 1.3  # 允许模拟冰源体积 / GM 体积 的上限

# V5优化时间划分
WARMUP_START = '2006-01-01'
WARMUP_END = '2008-12-31'
CALIB_START = '2009-01-01'
CALIB_END = '2017-12-31'
VALID_START = '2018-01-01'
VALID_END = '2020-12-31'
SIM_START = CALIB_START
SIM_END = VALID_END
OBS_MODE_OVERRIDE = 'auto'  # auto | full_year | may_oct

# =============================================================================
# 全局预加载数据
# =============================================================================
PREC_3D = None
TEMP_3D = None
ET_3D = None
LL_TEMP_3D = None
FLOW_ACC = None
PX_AREA = None
VALID_CELLS = None
Q_OBS_FULL = None
Q_OBS_OBJ_FULL = None  # 仅用于目标函数/指标评估的观测序列（例如冬季0值掩膜后）
Q_OBS_CALIB = None
Q_OBS_VALID = None
GLACIER_MELT_SERIES = None
GLACIER_MASK = None
GM_SERIES_SIM = None  # 与 SIM_DATES 对齐的 GM 融水序列（m3/s）
WARMUP_DAYS = None
SIM_DATES = None
SIM_MONTH = None
MAYJUN_MASK = None
MELT_MASK = None
WINTER_MASK = None
OBS_OBJ_MASK = None  # True where q_obs_obj is NaN (excluded from objective/metrics)
OBS_EVAL_MASK = None  # True where obs is used for objective/metrics (e.g., May-Oct for yearbook missing-season series)
OBS_EVAL_DESC = None
YEARBOOK_MODE = False
YEARBOOK_MAYOCT_FRAC = None
# 年鉴推算期的“月尺度”约束所需的预计算索引/统计量（与 SIM_DATES 对齐，均为 warmup-trim 后长度）
YEAR_VALUES = None        # 年列表（与 YEAR_SLICES 对应）
YEAR_SLICES = None        # [(start, end), ...]（end 为开区间）
OBS_YEAR_TOTAL = None     # 每年观测年总量（sum(q_obs)）
OBS_MONTH_SUM = None      # 每年每月观测月总量（shape: n_years x 12）
CALIB_MASK = None
VALID_MASK = None
AREA_COEF = None
ZONE_MID = None
ZONE_HIGH = None
CELL_SCALE = None  # 预计算的像元流量转换系数 (mm/d → m³/s)


# Logging
RUN_ID = None
LOG_FILE = None
PROGRESS_FILE = None
PROGRESS_FILE_MAIN = None
PROGRESS_FILE_REFINE = None
eval_count = 0
gen_count = 0  # differential_evolution 迭代代数（用于进度记录）
best_score = -np.inf
best_params = None
start_time = None
args = None

# Shared data for multiprocessing
SHARED_DATA = {}

# =============================================================================
# 参数边界（V6分区版 - 16参数）
# =============================================================================
# 说明：
# - ET 输入为 PM(PET)，ALPHA/ECORR 固定为 0。
# - 率定时加入必要的物理约束（如 K > K1 > K2；Muskingum K_MUSK 保持在系数非负域）。
param_names = [
    'TT', 'FC', 'BETA', 'LP',
    'RFCF', 'SFCF',
    'CFR', 'CWH',
    'CFMAX_mid', 'CFMAX_high',
    'K', 'K1', 'K2', 'PERC',
    'K_MUSK', 'ICE_FACTOR',
]

PARAM_BOUNDS = [
    (-2.0, 2.0),          # TT [degC]
    (100.0, 1500.0),      # FC [mm]
    (0.5, 4.0),           # BETA [-]
    (0.2, 1.0),           # LP [-]
    (0.8, 1.2),           # RFCF [-]
    (0.5, 1.5),           # SFCF [-]
    (0.0, 0.1),           # CFR [-]
    (0.01, 0.1),          # CWH [-]
    (2.0, 5.0),           # CFMAX_mid [mm/degC/day]
    (3.0, 8.0),           # CFMAX_high [mm/degC/day]
    (0.02, 0.5),          # K [1/day] (快速径流系数, 近似 K0)
    (0.01, 0.2),          # K1 [1/day] (中间响应/壤中流系数)
    (0.0005, 0.02),       # K2 [1/day] (基流系数; 需显著小于 K1)
    (0.01, 5.0),          # PERC [mm/day]
    (K_MUSK_MIN_SAFE, K_MUSK_MAX_SAFE),  # K_MUSK [day] (保证 Muskingum 系数非负)
    (1.2, 4.0),           # ICE_FACTOR [-] (裸冰度日因子倍率, >1)
]

FIXED = {
    'ALPHA': 0.0,   # ET 输入为 PM(PET)，不使用温度修正项
    'UZL': 30.0,    # 上层库快速出流阈值 [mm]（保持固定以降低参数不可辨识性）
    'X_MUSK': X_MUSK_FIXED,  # Muskingum X [-]
}

INIT_ST = np.array([0.0, 5.0, 0.0, 0.0, 0.0], dtype=np.float64)

# =============================================================================
# Numba-accelerated functions
# =============================================================================

@njit(cache=True, fastmath=True, nogil=True)
def hbv_cell(prec, temp, et, ll_temp, par, init_st, is_glacier, ice_factor):
    """Single cell HBV simulation with source-separated runoff."""
    n = len(prec)
    q_uz = np.zeros(n, dtype=np.float32)
    q_lz = np.zeros(n, dtype=np.float32)
    q_rain = np.zeros(n, dtype=np.float32)
    q_snow = np.zeros(n, dtype=np.float32)
    q_ice = np.zeros(n, dtype=np.float32)

    # Parse params
    tt, rfcf, sfcf, cfmax, cwh, cfr = par[0], par[1], par[2], par[3], par[4], par[5]
    fc, beta, e_corr, lp = par[6], par[7], par[8], par[9]
    k, k1, k2, uzl, perc = par[10], par[11], par[12], par[13], par[14]

    # Clip parameters to reasonable ranges
    fc = max(fc, 10.0)
    beta = max(min(beta, 10.0), 0.1)
    lp = max(min(lp, 0.99), 0.01)

    # Initial states
    sp, sm, uz, lz, wc = init_st[0], init_st[1], init_st[2], init_st[3], init_st[4]
    # Source tracking: start with no source attribution (will build up during warmup)
    uz_r, uz_s, uz_i = 0.0, 0.0, 0.0
    lz_r, lz_s, lz_i = 0.0, 0.0, 0.0

    for i in range(n):
        p = prec[i]
        t = temp[i]
        e = et[i]
        tm = ll_temp[i]

        # Skip invalid data
        if np.isnan(p) or np.isnan(t) or np.isnan(e):
            q_uz[i] = 0.0
            q_lz[i] = 0.0
            q_rain[i] = 0.0
            q_snow[i] = 0.0
            q_ice[i] = 0.0
            continue

        # Clip inputs to reasonable ranges
        p = max(p, 0.0)
        e = max(e, 0.0)

        # Precipitation
        if t <= tt:
            rf, sf = 0.0, p * sfcf
        else:
            rf, sf = p * rfcf, 0.0

        # Snow + glacier melt (two-stage within one day):
        # 1) snowmelt with CFMAX; 2) if glacier cell becomes (almost) snow-free, allow bare-ice melt with CFMAX*ICE_FACTOR.
        # This avoids "missing" melt energy on days when snow melts out and improves snow/ice source partition.
        melt = 0.0
        ice_melt = 0.0
        if t > tt:
            ddt = t - tt
            avail_snow = sp + sf
            melt_pot_snow = cfmax * ddt

            # Melt snow first (limited by available SWE)
            melt = min(melt_pot_snow, avail_snow)
            sp = max(avail_snow - melt, 0.0)

            # If glacier and snow melts out within the day, melt bare ice for the remaining time fraction.
            if is_glacier and (sp <= SWE_ICE_THRESHOLD_MM) and (melt_pot_snow > 0.0):
                f_snow = melt / (melt_pot_snow + EPS)  # 0..1
                if f_snow < 0.0:
                    f_snow = 0.0
                elif f_snow > 1.0:
                    f_snow = 1.0
                melt_pot_ice = (cfmax * ice_factor) * ddt
                ice_melt = max(melt_pot_ice * (1.0 - f_snow), 0.0)

            wc_int = wc + melt + rf + ice_melt
        else:
            refr = min(cfr * cfmax * (tt - t), wc + rf)
            sp = sp + sf + refr
            wc_int = max(wc - refr + rf, 0.0)

        # Clip snow pack
        sp = min(sp, 10000.0)

        if wc_int > cwh * sp:
            inf = wc_int - cwh * sp
            wc = cwh * sp
        else:
            inf = 0.0
            wc = wc_int

        # Source fractions (rain / snowmelt / icemelt)
        # 用 EPS 避免 0/0 在 fastmath 下触发 NaN；当 input_total=0 时三项本身均为 0。
        input_total = rf + melt + ice_melt
        den_in = input_total + EPS
        frac_r = rf / den_in
        frac_s = melt / den_in
        frac_i = ice_melt / den_in

        # Soil - with numerical stability
        sm_ratio = min(max(sm / fc, 0.0), 1.0)
        r = (sm_ratio ** beta) * inf
        # Recharge inherits source fractions from infiltration.
        r_r = r * frac_r
        r_s = r * frac_s
        r_i = r * frac_i

        ep_adj = max((1.0 + (t - tm) * e_corr) * e, 0.0)
        lp_fc = lp * fc
        if lp_fc > 0.001:
            ea = min(ep_adj, (sm / lp_fc) * ep_adj)
        else:
            ea = ep_adj
        ea = min(ea, sm)

        uz_r += r_r
        uz_s += r_s
        uz_i += r_i
        uz_int = uz_r + uz_s + uz_i
        sm = max(min(sm + inf - r - ea, fc), 0.0)

        # Response - with stability
        perc_actual = min(perc, uz_int)
        den_uz = uz_int + EPS
        frac_ur = uz_r / den_uz
        frac_us = uz_s / den_uz
        frac_ui = uz_i / den_uz
        uz_r -= perc_actual * frac_ur
        uz_s -= perc_actual * frac_us
        uz_i -= perc_actual * frac_ui
        lz_r += perc_actual * frac_ur
        lz_s += perc_actual * frac_us
        lz_i += perc_actual * frac_ui

        uz_int2 = uz_r + uz_s + uz_i

        q0 = k * max(uz_int2 - uzl, 0.0)
        q1 = k1 * uz_int2

        if q0 + q1 > uz_int2:
            q0 = uz_int2 * 0.67
            q1 = uz_int2 * 0.33
        den_uz2 = uz_int2 + EPS
        frac_ur2 = uz_r / den_uz2
        frac_us2 = uz_s / den_uz2
        frac_ui2 = uz_i / den_uz2

        q0_r = q0 * frac_ur2
        q0_s = q0 * frac_us2
        q0_i = q0 * frac_ui2
        q1_r = q1 * frac_ur2
        q1_s = q1 * frac_us2
        q1_i = q1 * frac_ui2

        uz_r = max(uz_r - (q0_r + q1_r), 0.0)
        uz_s = max(uz_s - (q0_s + q1_s), 0.0)
        uz_i = max(uz_i - (q0_i + q1_i), 0.0)

        lz_int = lz_r + lz_s + lz_i
        q2 = k2 * lz_int
        if q2 > lz_int:
            q2 = lz_int
        den_lz = lz_int + EPS
        frac_lr = lz_r / den_lz
        frac_ls = lz_s / den_lz
        frac_li = lz_i / den_lz

        q2_r = q2 * frac_lr
        q2_s = q2 * frac_ls
        q2_i = q2 * frac_li

        lz_r = max(lz_r - q2_r, 0.0)
        lz_s = max(lz_s - q2_s, 0.0)
        lz_i = max(lz_i - q2_i, 0.0)

        # Clip storages to prevent overflow
        uz_sum = uz_r + uz_s + uz_i
        if uz_sum > 10000.0:
            scale_uz = 10000.0 / uz_sum
            uz_r *= scale_uz
            uz_s *= scale_uz
            uz_i *= scale_uz

        lz_sum = lz_r + lz_s + lz_i
        if lz_sum > 10000.0:
            scale_lz = 10000.0 / lz_sum
            lz_r *= scale_lz
            lz_s *= scale_lz
            lz_i *= scale_lz

        q_uz[i] = q0 + q1
        q_lz[i] = q2
        q_rain[i] = q0_r + q1_r + q2_r
        q_snow[i] = q0_s + q1_s + q2_s
        q_ice[i] = q0_i + q1_i + q2_i

    return q_uz, q_lz, q_rain, q_snow, q_ice


@njit(cache=True, fastmath=True, nogil=True)
def run_all_cells(prec_3d, temp_3d, et_3d, ll_temp_3d, par_base, cfmax_grid, init_st,
                  valid_cells, cell_scale, glacier_mask, glacier_on, ice_factor):
    """Run HBV for all cells with zoned CFMAX (optimized with precomputed scale)."""
    n_cells = valid_cells.shape[0]
    ts = prec_3d.shape[2]

    # Aggregate discharge at outlet
    q_total = np.zeros(ts, dtype=np.float64)
    q_rain = np.zeros(ts, dtype=np.float64)
    q_snow = np.zeros(ts, dtype=np.float64)
    q_ice = np.zeros(ts, dtype=np.float64)

    for idx in range(n_cells):
        x = valid_cells[idx, 0]
        y = valid_cells[idx, 1]
        scale = cell_scale[idx]  # 预计算的转换系数

        prec = prec_3d[x, y, :]
        temp = temp_3d[x, y, :]
        et = et_3d[x, y, :]
        ll_temp = ll_temp_3d[x, y, :]

        # Get cell-specific CFMAX and create parameter array
        par = par_base.copy()
        par[3] = cfmax_grid[x, y]

        is_glacier = glacier_on and glacier_mask[x, y]
        q_uz, q_lz, q_r, q_s, q_i = hbv_cell(prec, temp, et, ll_temp, par, init_st, is_glacier, ice_factor)

        # Convert to m3/s and accumulate (使用预计算的scale)
        for t in range(ts):
            q_total[t] += (q_uz[t] + q_lz[t]) * scale
            q_rain[t] += q_r[t] * scale
            q_snow[t] += q_s[t] * scale
            q_ice[t] += q_i[t] * scale

    return q_total, q_rain, q_snow, q_ice


@njit(cache=True, fastmath=True, nogil=True)
def run_all_cells_glacier_snow(prec_3d, temp_3d, et_3d, ll_temp_3d, par_base, cfmax_grid, init_st,
                               valid_cells, cell_scale, glacier_mask, glacier_on, ice_factor):
    """Run HBV for all cells and additionally aggregate snowmelt from glacier cells only (pre-routing)."""
    n_cells = valid_cells.shape[0]
    ts = prec_3d.shape[2]

    q_total = np.zeros(ts, dtype=np.float64)
    q_rain = np.zeros(ts, dtype=np.float64)
    q_snow = np.zeros(ts, dtype=np.float64)
    q_ice = np.zeros(ts, dtype=np.float64)
    q_snow_glacier = np.zeros(ts, dtype=np.float64)

    for idx in range(n_cells):
        x = valid_cells[idx, 0]
        y = valid_cells[idx, 1]
        scale = cell_scale[idx]

        prec = prec_3d[x, y, :]
        temp = temp_3d[x, y, :]
        et = et_3d[x, y, :]
        ll_temp = ll_temp_3d[x, y, :]

        par = par_base.copy()
        par[3] = cfmax_grid[x, y]

        is_glacier = glacier_on and glacier_mask[x, y]
        q_uz, q_lz, q_r, q_s, q_i = hbv_cell(prec, temp, et, ll_temp, par, init_st, is_glacier, ice_factor)

        for t in range(ts):
            q_total[t] += (q_uz[t] + q_lz[t]) * scale
            q_rain[t] += q_r[t] * scale
            q_snow[t] += q_s[t] * scale
            q_ice[t] += q_i[t] * scale

        if is_glacier:
            for t in range(ts):
                q_snow_glacier[t] += q_s[t] * scale

    return q_total, q_rain, q_snow, q_ice, q_snow_glacier


@njit(cache=True, fastmath=True, nogil=True)
def muskingum_route(q_in, k, x):
    """Muskingum routing."""
    n = len(q_in)
    q_out = np.zeros(n, dtype=np.float64)
    q_out[0] = q_in[0]

    dt = 1.0
    denom = 2*k*(1-x) + dt
    c0 = (dt - 2*k*x) / denom
    c1 = (dt + 2*k*x) / denom
    c2 = (2*k*(1-x) - dt) / denom

    for i in range(1, n):
        q_out[i] = c0 * q_in[i] + c1 * q_in[i-1] + c2 * q_out[i-1]
        if q_out[i] < 0:
            q_out[i] = 0.0

    return q_out


@njit(cache=True, nogil=True)
def nse_numba(obs, sim):
    """Numba-optimized NSE."""
    n = len(obs)
    sum_obs = 0.0
    count = 0

    for i in range(n):
        if not (np.isnan(obs[i]) or np.isnan(sim[i])):
            sum_obs += obs[i]
            count += 1

    if count == 0:
        return -999.0

    mean_obs = sum_obs / count
    ss_err = 0.0
    ss_tot = 0.0

    for i in range(n):
        if not (np.isnan(obs[i]) or np.isnan(sim[i])):
            ss_err += (obs[i] - sim[i]) ** 2
            ss_tot += (obs[i] - mean_obs) ** 2

    if ss_tot == 0:
        return -999.0

    return 1.0 - ss_err / ss_tot


def warmup_jit():
    """Pre-compile all JIT functions."""
    print("  Warming up Numba JIT...")

    n = 10
    prec = np.random.rand(n).astype(np.float32)
    temp = np.random.rand(n).astype(np.float32) * 10
    et = np.random.rand(n).astype(np.float32)
    ll_temp = np.random.rand(n).astype(np.float32) * 5

    # TT, RFCF, SFCF, CFMAX, CWH, CFR, FC, BETA, ALPHA(ECORR), LP, K, K1, K2, UZL, PERC
    par = np.array([0.0, 1.0, 1.0, 3.0, 0.05, 0.05, 200.0, 1.0, 0.0, 0.5,
                    0.2, 0.05, 0.005, 30.0, 2.0], dtype=np.float64)
    init_st = np.array([0.0, 5.0, 5.0, 5.0, 0.0], dtype=np.float64)

    _ = hbv_cell(prec, temp, et, ll_temp, par, init_st, True, 1.5)

    rows, cols = 5, 5
    prec_3d = np.random.rand(rows, cols, n).astype(np.float32)
    temp_3d = np.random.rand(rows, cols, n).astype(np.float32) * 10
    et_3d = np.random.rand(rows, cols, n).astype(np.float32)
    ll_temp_3d = np.random.rand(rows, cols, n).astype(np.float32) * 5
    valid_cells = np.array([[i, j] for i in range(rows) for j in range(cols)], dtype=np.int64)
    cell_scale = np.ones(rows * cols, dtype=np.float64) / 86.4
    cfmax_grid = np.ones((rows, cols), dtype=np.float64) * 4.0

    glacier_mask = np.zeros((rows, cols), dtype=np.bool_)
    _ = run_all_cells(prec_3d, temp_3d, et_3d, ll_temp_3d, par, cfmax_grid, init_st,
                      valid_cells, cell_scale, glacier_mask, True, 1.5)
    _ = run_all_cells_glacier_snow(prec_3d, temp_3d, et_3d, ll_temp_3d, par, cfmax_grid, init_st,
                                   valid_cells, cell_scale, glacier_mask, True, 1.5)


    q = np.random.rand(n).astype(np.float64)
    _ = muskingum_route(q, 1.0, 0.2)

    _ = nse_numba(q, q)

    print("  JIT compilation complete!")


# =============================================================================
# Data loading
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="HBV calibration V6 with zoned CFMAX")
    parser.add_argument("--maxiter", type=int, default=80)
    parser.add_argument("--popsize", type=int, default=12)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--log-every", type=int, default=10)
    # SciPy DE 的收敛判据依赖种群能量方差：tol 过大会在 1-2 代就“看似收敛”提前停止。
    # 这里默认更严格一些，避免 progress.csv/progress_refine.csv 只有 1-2 行。
    parser.add_argument("--tol", type=float, default=0.001)
    # 默认关闭 polish：SciPy 的局部优化阶段可能长时间无输出，易被误认为“卡死”
    parser.add_argument("--polish", action="store_true", default=False)
    parser.add_argument("--no-polish", dest="polish", action="store_false")
    parser.add_argument("--prec-source", choices=["mswep", "cmfd"], default="mswep")
    parser.add_argument("--prec-dir", type=str, default=None, help="Override precipitation directory")
    parser.add_argument("--glacier-mode", choices=["inline", "series", "off"], default="inline")
    parser.add_argument("--fill-nan", dest="fill_nan", action="store_true", default=True)
    parser.add_argument("--no-fill-nan", dest="fill_nan", action="store_false")
    parser.add_argument("--quick-test", action="store_true", default=False)
    parser.add_argument("--quick-days", type=int, default=30)
    parser.add_argument("--analyze-run", type=str, default=None, help="Analyze an existing run directory (no calibration).")
    # 默认：启用 GM 对照（仅用于诊断/独立验证，不参与优化）；可用 --skip-gm 关闭。
    parser.add_argument("--skip-gm", dest="skip_gm", action="store_true", default=None,
                        help="Skip reading GM_*.tif glacier melt series in diagnostics.")
    parser.add_argument("--use-gm", dest="skip_gm", action="store_false",
                        help="Enable reading GM_*.tif glacier melt series in diagnostics.")
    return parser.parse_args()


def setup_logging():
    global RUN_ID, LOG_FILE, PROGRESS_FILE, PROGRESS_FILE_MAIN, start_time, gen_count
    RUN_ID = datetime.now().strftime("%Y%m%d_%H%M%S")
    os.makedirs(LOG_DIR, exist_ok=True)
    LOG_FILE = os.path.join(LOG_DIR, f"v6_{RUN_ID}.log")
    PROGRESS_FILE = os.path.join(LOG_DIR, f"v6_progress_{RUN_ID}.csv")
    PROGRESS_FILE_MAIN = PROGRESS_FILE
    start_time = time.time()
    gen_count = 0
    with open(PROGRESS_FILE, 'w') as f:
        f.write("timestamp,gen,nse_cal,nse_log_cal,nse_val,nse_log_val,pbias_cal,pbias_val,obj,convergence,elapsed_sec\n")


def log_msg(msg):
    print(msg)
    if LOG_FILE:
        with open(LOG_FILE, 'a', encoding='utf-8') as f:
            f.write(msg + "\n")


def append_progress(gen, nse_cal, nse_log_cal, nse_val, nse_log_val, pbias_cal, pbias_val, obj, convergence):
    elapsed = time.time() - start_time
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(PROGRESS_FILE, 'a') as f:
        f.write(
            f"{ts},{gen},{nse_cal:.6f},{nse_log_cal:.6f},{nse_val:.6f},{nse_log_val:.6f},"
            f"{pbias_cal:.3f},{pbias_val:.3f},{obj:.6f},{convergence:.6f},{elapsed:.1f}\n"
        )


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
    """Load time series of GeoTIFFs."""
    dates = pd.date_range(start_date, end_date, freq='D')

    all_files = glob(os.path.join(directory, '*.tif'))
    if not all_files:
        raise ValueError(f"No .tif files in {directory}")

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


def load_glacier_melt_series(start_date, end_date):
    # Cache: avoid re-reading thousands of daily rasters across repeated runs/diagnostics.
    try:
        os.makedirs(CACHE_DIR, exist_ok=True)
        cache_key = f"gm_series_{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}_f{GLACIER_MELT_FACTOR:.3f}.npy"
        cache_file = os.path.join(CACHE_DIR, cache_key)
        if os.path.exists(cache_file):
            arr = np.load(cache_file)
            return arr.astype(np.float64, copy=False)
    except Exception:
        cache_file = None

    melt_files = sorted(glob(os.path.join(GLACIER_MELT_DIR, 'GM_*.tif')))
    if not melt_files:
        return None

    records = {}
    row_areas = None

    for melt_file in melt_files:
        name = os.path.basename(melt_file)
        try:
            date_str = name.replace('GM_', '').replace('.tif', '')
            date = datetime.strptime(date_str, '%Y.%m.%d').date()
        except:
            continue

        if date < start_date or date > end_date:
            continue

        with rasterio.open(melt_file) as src:
            melt = src.read(1)
            nodata = src.nodata
            if row_areas is None:
                row_areas = compute_row_areas_km2(src.transform, src.height)

        if nodata is not None:
            melt = np.where(melt == nodata, 0, melt)
        melt = np.where(melt > 0, melt, 0)
        melt_m3s = float((melt * row_areas[:, None]).sum() / 86.4) * GLACIER_MELT_FACTOR
        records[date] = melt_m3s

    if not records:
        return None

    idx = pd.date_range(start_date, end_date, freq='D')
    series = np.array([records.get(d.date(), 0.0) for d in idx], dtype=np.float64)
    if cache_file:
        try:
            np.save(cache_file, series)
        except Exception:
            pass
    return series


def detect_obs_eval_mask(q_obs_full):
    full_mask = np.ones(len(SIM_DATES), dtype=bool)
    may_oct_mask = (SIM_MONTH >= 5) & (SIM_MONTH <= 10)
    mode = str(OBS_MODE_OVERRIDE or 'auto').strip().lower()

    if mode == 'full_year':
        return full_mask, 'full_year_forced', False, None
    if mode == 'may_oct':
        return may_oct_mask, 'mayoct_forced', True, YEARBOOK_MAYOCT_FRAC_DEFAULT
    if mode != 'auto':
        print(f"      [WARN] Unknown OBS_MODE_OVERRIDE={OBS_MODE_OVERRIDE!r}, fallback to auto.")

    try:
        df_frac = pd.DataFrame({'date': SIM_DATES, 'q': q_obs_full})
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
            return may_oct_mask, f"mayoct_only_yearbook(frac={frac_mean:.3f}±{frac_std:.3g})", True, float(frac_mean)
    except Exception as exc:
        print(f"      [WARN] Obs mode detection failed, fallback to full-year eval: {exc}")

    return full_mask, 'full_year', False, None


def load_all_data(end_date_override=None, skip_obs=False):
    global PREC_3D, TEMP_3D, ET_3D, LL_TEMP_3D, FLOW_ACC, PX_AREA
    global VALID_CELLS, Q_OBS_FULL, Q_OBS_OBJ_FULL, Q_OBS_CALIB, Q_OBS_VALID, GLACIER_MELT_SERIES, GLACIER_MASK, WARMUP_DAYS, AREA_COEF
    global SIM_DATES, CALIB_MASK, VALID_MASK
    global SIM_MONTH, MAYJUN_MASK, MELT_MASK, WINTER_MASK
    global OBS_OBJ_MASK, OBS_EVAL_MASK, OBS_EVAL_DESC
    global YEARBOOK_MODE, YEARBOOK_MAYOCT_FRAC
    global YEAR_VALUES, YEAR_SLICES, OBS_YEAR_TOTAL, OBS_MONTH_SUM
    global ZONE_MID, ZONE_HIGH

    load_end = SIM_END
    if end_date_override:
        load_end = end_date_override

    print("="*70)
    print("Loading data...")
    print("="*70)

    print("[1/5] Loading precipitation...")
    PREC_3D, transform = load_raster_stack(PREC_DIR, WARMUP_START, load_end)
    rows, cols, ts = PREC_3D.shape
    print(f"      Shape: {rows} x {cols} x {ts}")

    print("[2/5] Loading temperature...")
    TEMP_3D, _ = load_raster_stack(TEMP_DIR, WARMUP_START, load_end)

    print("[3/5] Loading ET...")
    ET_3D, _ = load_raster_stack(EVAP_DIR, WARMUP_START, load_end)

    # Long-term mean temperature
    LL_TEMP_3D = np.nanmean(TEMP_3D, axis=2, keepdims=True)
    LL_TEMP_3D = np.broadcast_to(LL_TEMP_3D, TEMP_3D.shape).astype(np.float32)

    # Fill NaNs to avoid empty cells producing zero flow
    if args.fill_nan:
        PREC_3D = np.nan_to_num(PREC_3D, nan=0.0)
        ET_3D = np.nan_to_num(ET_3D, nan=0.0)
        TEMP_3D = np.where(np.isnan(TEMP_3D), LL_TEMP_3D, TEMP_3D)

    print("[4/5] Loading flow accumulation...")
    with rasterio.open(FLOW_ACC_PATH) as src:
        FLOW_ACC = src.read(1).astype(np.float32)
        nodata = src.nodata
        if nodata is not None:
            FLOW_ACC[FLOW_ACC == nodata] = np.nan

    valid_mask = ~np.isnan(FLOW_ACC)
    VALID_CELLS = np.argwhere(valid_mask).astype(np.int64)
    print(f"      Valid cells: {len(VALID_CELLS)}")

    # Load elevation zones (V6)
    print("[4.5/5] Loading elevation zones...")
    zone_mid_file = os.path.join(PROJECT_ROOT, 'data/gis/elevation_zone_mid.tif')
    zone_high_file = os.path.join(PROJECT_ROOT, 'data/gis/elevation_zone_high.tif')
    with rasterio.open(zone_mid_file) as src:
        ZONE_MID = src.read(1).astype(bool)
    with rasterio.open(zone_high_file) as src:
        ZONE_HIGH = src.read(1).astype(bool)
    print(f"      Mid zone: {ZONE_MID.sum()} cells, High zone: {ZONE_HIGH.sum()} cells")

    # Pixel areas
    PX_AREA = np.zeros((rows, cols), dtype=np.float64)
    row_areas = compute_row_areas_km2(transform, rows)
    for i in range(rows):
        PX_AREA[i, :] = row_areas[i]

    px_tot_area = np.sum(PX_AREA[valid_mask])
    AREA_COEF = CATCHMENT_AREA / px_tot_area

    # Precompute cell scale factors for speed optimization (mm/d → m³/s)
    global CELL_SCALE
    CELL_SCALE = (PX_AREA[VALID_CELLS[:,0], VALID_CELLS[:,1]] * AREA_COEF / 86.4).astype(np.float64)
    print(f"      Cell scale factors computed: {len(CELL_SCALE)} cells")

    print("[5/5] Loading observed discharge...")
    
    # Reset yearbook-mode flags (safe for repeated runs within one Python session)
    YEARBOOK_MODE = False
    YEARBOOK_MAYOCT_FRAC = None

    # SIM_DATES matches the TRIMMED simulation output (after removing warmup)
    # run_simulation returns q_sim[WARMUP_DAYS:], so SIM_DATES should start from SIM_START (=CALIB_START)
    sim_end_date = load_end if load_end else SIM_END
    SIM_DATES = pd.date_range(start=SIM_START, end=sim_end_date, freq='D')

    # 预计算月份掩膜（用于目标函数的季节性约束）
    SIM_MONTH = SIM_DATES.month.values.astype(np.int16)
    MAYJUN_MASK = (SIM_MONTH == 5) | (SIM_MONTH == 6)
    MELT_MASK = (SIM_MONTH >= 4) & (SIM_MONTH <= 9)
    WINTER_MASK = (SIM_MONTH == 12) | (SIM_MONTH <= 2)
    
    if not skip_obs:
        if not os.path.exists(OBS_FILE):
             print(f"      [WARN] Obs file not found: {OBS_FILE}")
             Q_OBS_FULL = np.zeros(len(SIM_DATES))
             Q_OBS_OBJ_FULL = Q_OBS_FULL.copy()
             OBS_EVAL_MASK = np.ones(len(SIM_DATES), dtype=bool)
             OBS_EVAL_DESC = "no_obs_file"
             OBS_OBJ_MASK = np.zeros(len(SIM_DATES), dtype=bool)
        else:
            df_obs = pd.read_csv(OBS_FILE)
            df_obs['date'] = pd.to_datetime(df_obs['date'])
            df_obs.set_index('date', inplace=True)
            Q_OBS_FULL = df_obs.reindex(SIM_DATES)['discharge (m3/s)'].values.astype(np.float64)

            OBS_EVAL_MASK, OBS_EVAL_DESC, YEARBOOK_MODE, YEARBOOK_MAYOCT_FRAC = detect_obs_eval_mask(Q_OBS_FULL)
            if YEARBOOK_MODE:
                print(f"      Obs mode detected: {OBS_EVAL_DESC}.")
                print("      Objective/metrics use observed window only: May-Oct (5-10).")
                target_frac = YEARBOOK_MAYOCT_FRAC_DEFAULT if YEARBOOK_MAYOCT_FRAC is None else float(YEARBOOK_MAYOCT_FRAC)
                min_frac = max(0.0, target_frac - YEARBOOK_MAYOCT_FRAC_SLACK)
                print(f"      Seasonality prior enabled: require simulated May-Oct fraction >= {min_frac:.3f} (soft penalty).")
            else:
                print(f"      Obs mode detected: {OBS_EVAL_DESC}.")
                print("      Objective/metrics use full-year daily observations.")

            Q_OBS_OBJ_FULL = Q_OBS_FULL.copy()
            OBS_OBJ_MASK = (~OBS_EVAL_MASK).astype(bool)
            if np.any(OBS_OBJ_MASK):
                Q_OBS_OBJ_FULL[OBS_OBJ_MASK] = np.nan

        CALIB_MASK = (SIM_DATES >= pd.to_datetime(CALIB_START)) & (SIM_DATES <= pd.to_datetime(CALIB_END))
        VALID_MASK = (SIM_DATES >= pd.to_datetime(VALID_START)) & (SIM_DATES <= pd.to_datetime(VALID_END))
        if Q_OBS_OBJ_FULL is None:
            Q_OBS_OBJ_FULL = Q_OBS_FULL
        Q_OBS_CALIB = Q_OBS_OBJ_FULL[CALIB_MASK]
        Q_OBS_VALID = Q_OBS_OBJ_FULL[VALID_MASK]
        print(f"      Observed mean (calib): {np.nanmean(Q_OBS_CALIB):.2f} m3/s")
    else:
        # Dummy obs for quick test
        Q_OBS_FULL = np.zeros(len(SIM_DATES))
        Q_OBS_OBJ_FULL = Q_OBS_FULL.copy()
        OBS_EVAL_MASK = np.ones(len(SIM_DATES), dtype=bool)
        OBS_EVAL_DESC = "skip_obs"
        OBS_OBJ_MASK = np.zeros(len(SIM_DATES), dtype=bool)
        Q_OBS_CALIB = np.array([])
        Q_OBS_VALID = np.array([])
        CALIB_MASK = np.zeros(len(SIM_DATES), dtype=bool)
        VALID_MASK = np.zeros(len(SIM_DATES), dtype=bool)

    # 预计算：逐年索引与观测月尺度统计（用于年鉴推算期的“月体积”约束项）
    # 注意：这些统计与 SIM_DATES 对齐（已经是 warmup-trim 后的模拟窗口）。
    try:
        years_arr = SIM_DATES.year.values.astype(np.int32)
        YEAR_VALUES = []
        YEAR_SLICES = []
        if len(years_arr) > 0:
            s0 = 0
            for i in range(1, len(years_arr)):
                if years_arr[i] != years_arr[i - 1]:
                    YEAR_VALUES.append(int(years_arr[s0]))
                    YEAR_SLICES.append((int(s0), int(i)))
                    s0 = i
            YEAR_VALUES.append(int(years_arr[s0]))
            YEAR_SLICES.append((int(s0), int(len(years_arr))))
        OBS_YEAR_TOTAL = None
        OBS_MONTH_SUM = None
        if Q_OBS_FULL is not None and len(Q_OBS_FULL) == len(SIM_DATES):
            OBS_YEAR_TOTAL = np.zeros(len(YEAR_VALUES), dtype=np.float64)
            OBS_MONTH_SUM = np.zeros((len(YEAR_VALUES), 12), dtype=np.float64)
            for yi, (s, e) in enumerate(YEAR_SLICES):
                qy = np.asarray(Q_OBS_FULL[s:e], dtype=np.float64)
                my = np.asarray(SIM_MONTH[s:e], dtype=np.int16)
                OBS_YEAR_TOTAL[yi] = float(np.nansum(qy))
                for m in range(1, 13):
                    OBS_MONTH_SUM[yi, m - 1] = float(np.nansum(qy[my == m]))
    except Exception as e:
        print(f"      [WARN] Yearbook monthly precompute failed (skip monthly constraints): {e}")
        YEAR_VALUES = None
        YEAR_SLICES = None
        OBS_YEAR_TOTAL = None
        OBS_MONTH_SUM = None

    # Glacier mask
    if os.path.exists(GLACIER_MASK_PATH):
        with rasterio.open(GLACIER_MASK_PATH) as src:
            glacier_raw = src.read(1)
            nodata = src.nodata
            if nodata is not None:
                glacier_raw = np.where(glacier_raw == nodata, 0, glacier_raw)
        GLACIER_MASK = glacier_raw.astype(bool)
        print(f"      Glacier cells: {GLACIER_MASK.sum()}")
    else:
        GLACIER_MASK = np.zeros((rows, cols), dtype=bool)
        print("      Glacier mask missing, using empty mask")

    WARMUP_DAYS = (pd.to_datetime(CALIB_START) - pd.to_datetime(WARMUP_START)).days
    print(f"      Warmup days: {WARMUP_DAYS}")

    # Glacier melt series (GM_*.tif)
    # - series 模式：作为外加冰川融水直接叠加到出口流量
    # - inline/off 模式：默认不参与模拟，但可用于目标函数/诊断的独立一致性约束
    global GM_SERIES_SIM
    GM_SERIES_SIM = None
    GLACIER_MELT_SERIES = np.zeros(ts, dtype=np.float64)
    if os.path.isdir(GLACIER_MELT_DIR) and (args.glacier_mode == "series" or (OBJ_W_GM_VOL > 0) or (args.skip_gm is False)):
        start_date = datetime.strptime(WARMUP_START, '%Y-%m-%d').date()
        end_date = datetime.strptime(load_end, '%Y-%m-%d').date()
        gm_full = load_glacier_melt_series(start_date, end_date)
        if gm_full is None:
            gm_full = np.zeros(ts, dtype=np.float64)
        else:
            gm_full = gm_full.astype(np.float64, copy=False)

        # 与 SIM_DATES 对齐（去除 warmup）
        if WARMUP_DAYS < len(gm_full):
            GM_SERIES_SIM = gm_full[WARMUP_DAYS:WARMUP_DAYS + len(SIM_DATES)]

        if args.glacier_mode == "series":
            GLACIER_MELT_SERIES = gm_full
            print(f"      Glacier melt mean (GM): {np.mean(GLACIER_MELT_SERIES):.2f} m3/s")

    # Store all data in SHARED_DATA for multiprocessing
    global SHARED_DATA
    SHARED_DATA = {
        'PREC_3D': PREC_3D,
        'TEMP_3D': TEMP_3D,
        'ET_3D': ET_3D,
        'LL_TEMP_3D': LL_TEMP_3D,
        'VALID_CELLS': VALID_CELLS,
        'CELL_SCALE': CELL_SCALE,
        'GLACIER_MASK': GLACIER_MASK,
        'GLACIER_MELT_SERIES': GLACIER_MELT_SERIES,
        'GM_SERIES_SIM': GM_SERIES_SIM,
        'ZONE_MID': ZONE_MID,
        'ZONE_HIGH': ZONE_HIGH,
        'Q_OBS_FULL': Q_OBS_FULL,
        'Q_OBS_CALIB': Q_OBS_CALIB,
        'Q_OBS_VALID': Q_OBS_VALID,
        'CALIB_MASK': CALIB_MASK,
        'VALID_MASK': VALID_MASK,
        'WARMUP_DAYS': WARMUP_DAYS,
        'SIM_DATES': SIM_DATES,
        'INIT_ST': INIT_ST,
        'FIXED': FIXED,
        'args_glacier_mode': args.glacier_mode,
    }

    print("="*70)


def run_simulation(opt_params):
    """Run HBV simulation with optimized params (V6 with zoned CFMAX + glacier mode)."""
    TT, FC, BETA, LP, RFCF, SFCF, CFR, CWH, CFMAX_mid, CFMAX_high, K, K1, K2, PERC, K_MUSK, ICE_FACTOR = opt_params

    # Create CFMAX grid based on elevation zones
    cfmax_grid = np.zeros((PREC_3D.shape[0], PREC_3D.shape[1]), dtype=np.float64)
    cfmax_grid[ZONE_MID] = CFMAX_mid
    cfmax_grid[ZONE_HIGH] = CFMAX_high
    # Fallback for any cells not covered by the zone rasters (use mid-elevation CFMAX).
    cfmax_grid[~(ZONE_MID | ZONE_HIGH)] = CFMAX_mid

    # Build base parameter array (15 params for hbv_cell, CFMAX will be set per cell)
    par_base = np.array([
        TT, RFCF, SFCF, 0.0, CWH, CFR,
        FC, BETA, FIXED['ALPHA'], LP, K, K1, K2, FIXED['UZL'], PERC
    ], dtype=np.float64)

    glacier_on = args.glacier_mode == "inline"

    q_total, q_rain, q_snow, q_ice = run_all_cells(
        PREC_3D, TEMP_3D, ET_3D, LL_TEMP_3D,
        par_base, cfmax_grid, INIT_ST, VALID_CELLS, CELL_SCALE,
        GLACIER_MASK, glacier_on, ICE_FACTOR
    )

    # Apply Muskingum routing for catchment response
    q_routed = muskingum_route(q_total, K_MUSK, FIXED['X_MUSK'])
    q_routed_rain = muskingum_route(q_rain, K_MUSK, FIXED['X_MUSK'])
    q_routed_snow = muskingum_route(q_snow, K_MUSK, FIXED['X_MUSK'])
    q_routed_ice = muskingum_route(q_ice, K_MUSK, FIXED['X_MUSK'])

    # Glacier handling
    if args.glacier_mode == "series":
        min_len = min(len(q_routed), len(GLACIER_MELT_SERIES))
        q_final = q_routed[:min_len] + GLACIER_MELT_SERIES[:min_len]
        q_rain = q_routed_rain[:min_len]
        q_snow = q_routed_snow[:min_len]
        q_ice = GLACIER_MELT_SERIES[:min_len]
    elif args.glacier_mode == "off":
        q_final = q_routed
        q_rain = q_routed_rain
        q_snow = q_routed_snow
        q_ice = np.zeros_like(q_final)
    else:
        q_final = q_routed
        q_rain = q_routed_rain
        q_snow = q_routed_snow
        q_ice = q_routed_ice

    # Return calibration period
    return (
        q_final[WARMUP_DAYS:],
        q_rain[WARMUP_DAYS:],
        q_snow[WARMUP_DAYS:],
        q_ice[WARMUP_DAYS:],
    )


def run_simulation_detailed(opt_params):
    """Forward run with extra cryosphere diagnostics (computed only once, not used in objective)."""
    TT, FC, BETA, LP, RFCF, SFCF, CFR, CWH, CFMAX_mid, CFMAX_high, K, K1, K2, PERC, K_MUSK, ICE_FACTOR = opt_params

    cfmax_grid = np.zeros((PREC_3D.shape[0], PREC_3D.shape[1]), dtype=np.float64)
    cfmax_grid[ZONE_MID] = CFMAX_mid
    cfmax_grid[ZONE_HIGH] = CFMAX_high
    cfmax_grid[~(ZONE_MID | ZONE_HIGH)] = CFMAX_mid

    par_base = np.array([
        TT, RFCF, SFCF, 0.0, CWH, CFR,
        FC, BETA, FIXED['ALPHA'], LP, K, K1, K2, FIXED['UZL'], PERC
    ], dtype=np.float64)

    glacier_on = args.glacier_mode == "inline"

    q_total_raw, q_rain_raw, q_snow_raw, q_ice_raw, q_snow_glacier_raw = run_all_cells_glacier_snow(
        PREC_3D, TEMP_3D, ET_3D, LL_TEMP_3D,
        par_base, cfmax_grid, INIT_ST, VALID_CELLS, CELL_SCALE,
        GLACIER_MASK, glacier_on, ICE_FACTOR
    )

    # Routing (linear, so we can safely route sub-components and keep closure)
    q_total = muskingum_route(q_total_raw, K_MUSK, FIXED['X_MUSK'])
    q_rain = muskingum_route(q_rain_raw, K_MUSK, FIXED['X_MUSK'])
    q_snow = muskingum_route(q_snow_raw, K_MUSK, FIXED['X_MUSK'])
    q_ice = muskingum_route(q_ice_raw, K_MUSK, FIXED['X_MUSK'])
    q_snow_glacier = muskingum_route(q_snow_glacier_raw, K_MUSK, FIXED['X_MUSK'])
    q_glacier_total = q_snow_glacier + q_ice

    # Glacier handling
    if args.glacier_mode == "series":
        # 外部冰川融水序列通常代表“冰川产流总量”，此时不再对 q_ice/q_glacier_total 做 inline 分解。
        min_len = min(len(q_total), len(GLACIER_MELT_SERIES))
        q_final = q_total[:min_len] + GLACIER_MELT_SERIES[:min_len]
        q_rain = q_rain[:min_len]
        q_snow = q_snow[:min_len]
        q_ice = GLACIER_MELT_SERIES[:min_len]
        q_snow_glacier = np.zeros_like(q_ice)
        q_glacier_total = q_ice.copy()
        q_total_raw = q_total_raw[:min_len]
        q_rain_raw = q_rain_raw[:min_len]
        q_snow_raw = q_snow_raw[:min_len]
        q_ice_raw = GLACIER_MELT_SERIES[:min_len].copy()
        q_snow_glacier_raw = np.zeros_like(q_ice_raw)
    elif args.glacier_mode == "off":
        q_final = q_total
        q_ice = np.zeros_like(q_final)
        q_snow_glacier = np.zeros_like(q_final)
        q_glacier_total = np.zeros_like(q_final)
    else:
        q_final = q_total

    # Align with SIM_DATES (simulation outputs always drop warmup)
    return (
        q_final[WARMUP_DAYS:],
        q_rain[WARMUP_DAYS:],
        q_snow[WARMUP_DAYS:],
        q_ice[WARMUP_DAYS:],
        q_snow_glacier[WARMUP_DAYS:],
        q_glacier_total[WARMUP_DAYS:],
        # raw (pre-routing) for diagnostics/comparison with GM (also warmup-trimmed)
        q_total_raw[WARMUP_DAYS:],
        q_rain_raw[WARMUP_DAYS:],
        q_snow_raw[WARMUP_DAYS:],
        q_ice_raw[WARMUP_DAYS:],
        q_snow_glacier_raw[WARMUP_DAYS:],
    )


def _init_worker(shared_data):
    """Initialize worker process with shared data."""
    global PREC_3D, TEMP_3D, ET_3D, LL_TEMP_3D, VALID_CELLS, CELL_SCALE
    global GLACIER_MASK, GLACIER_MELT_SERIES, GM_SERIES_SIM, ZONE_MID, ZONE_HIGH
    global Q_OBS_CALIB, CALIB_MASK, WARMUP_DAYS, INIT_ST, FIXED, args
    
    PREC_3D = shared_data['PREC_3D']
    TEMP_3D = shared_data['TEMP_3D']
    ET_3D = shared_data['ET_3D']
    LL_TEMP_3D = shared_data['LL_TEMP_3D']
    VALID_CELLS = shared_data['VALID_CELLS']
    CELL_SCALE = shared_data['CELL_SCALE']
    GLACIER_MASK = shared_data['GLACIER_MASK']
    GLACIER_MELT_SERIES = shared_data['GLACIER_MELT_SERIES']
    GM_SERIES_SIM = shared_data.get('GM_SERIES_SIM', None)
    ZONE_MID = shared_data['ZONE_MID']
    ZONE_HIGH = shared_data['ZONE_HIGH']
    Q_OBS_CALIB = shared_data['Q_OBS_CALIB']
    CALIB_MASK = shared_data['CALIB_MASK']
    WARMUP_DAYS = shared_data['WARMUP_DAYS']
    INIT_ST = shared_data['INIT_ST']
    FIXED = shared_data['FIXED']
    
    # Create a simple args-like object for glacier_mode
    class Args:
        pass
    args = Args()
    args.glacier_mode = shared_data['args_glacier_mode']


def _seasonality_penalty(q_sim):
    """
    年鉴缺测/推算期模式下的弱季节性约束。

    - 目标函数主拟合仅使用 5-10 月实测窗口，避免被推算期的“常数日值”抬高 NSE。
    - 但仍需避免参数把水量不合理地转移到冬季（11-4 月），因此用 5-10 月占比做软约束。
    """
    if not YEARBOOK_MODE:
        return 0.0, None
    if OBJ_W_SEASON <= 0:
        return 0.0, None
    if OBS_EVAL_MASK is None:
        return 0.0, None

    q = np.asarray(q_sim, dtype=np.float64)
    tot = float(np.nansum(q))
    if (not np.isfinite(tot)) or tot <= 0.0:
        return 0.0, None

    frac_mayoct = float(np.nansum(q[OBS_EVAL_MASK]) / (tot + EPS))
    target = YEARBOOK_MAYOCT_FRAC_DEFAULT
    if YEARBOOK_MAYOCT_FRAC is not None and np.isfinite(YEARBOOK_MAYOCT_FRAC):
        target = float(YEARBOOK_MAYOCT_FRAC)

    min_frac = max(0.0, target - YEARBOOK_MAYOCT_FRAC_SLACK)
    if frac_mayoct >= min_frac:
        return 0.0, frac_mayoct

    pen = float(OBJ_W_SEASON * (min_frac - frac_mayoct) ** 2)
    return pen, frac_mayoct


def _offseason_monthly_penalty(q_sim):
    """
    年鉴推算期（11-4 月）按月体积约束：
    - 不把“推算日值(常数/0)”当作真实逐日过程去拟合 NSE；
    - 但充分利用其提供的“月尺度体积/季节分配”信息，抑制冬季高流量并约束 4 月回升。
    """
    if not YEARBOOK_MODE:
        return 0.0, None
    if OBJ_W_OFFSEASON_MONTHLY <= 0:
        return 0.0, None
    if (YEAR_SLICES is None) or (OBS_YEAR_TOTAL is None) or (OBS_MONTH_SUM is None) or (SIM_MONTH is None):
        return 0.0, None

    q = np.asarray(q_sim, dtype=np.float64)
    pen = 0.0
    worst_excess_zero = 0.0  # 目标为 0 的月份（12-2）超额占比（每月）

    for yi, (s, e) in enumerate(YEAR_SLICES):
        obs_year = float(OBS_YEAR_TOTAL[yi]) if OBS_YEAR_TOTAL is not None else 0.0
        if (not np.isfinite(obs_year)) or obs_year <= 0.0:
            continue

        qy = q[s:e]
        my = np.asarray(SIM_MONTH[s:e], dtype=np.int16)

        for m in OFFSEASON_MONTHS:
            sim_m = float(np.nansum(qy[my == m]))
            obs_m = float(OBS_MONTH_SUM[yi, m - 1])

            # 对目标为 0 的月份：允许少量基流（占全年体积的阈值），超过部分才惩罚
            if obs_m <= 0.0:
                sim_frac = sim_m / (obs_year + EPS)
                excess = sim_frac - OFFSEASON_ZERO_FRAC_SLACK
                if excess > 0.0:
                    # 目标为 0 的月份对物理约束更关键，权重略高
                    w = 4.0
                    pen += w * excess * excess
                    if excess > worst_excess_zero:
                        worst_excess_zero = excess
                continue

            # 对目标为正的月份：按“占全年体积的相对误差”惩罚（避免单位依赖）
            rel = (sim_m - obs_m) / (obs_year + EPS)
            # 4 月在年鉴推算中权重更大（回升段），这里给予略高权重，提升 4 月/5 月衔接
            w = 2.0 if m == 4 else 1.0
            pen += w * rel * rel

    return float(OBJ_W_OFFSEASON_MONTHLY * pen), float(worst_excess_zero)


def _winter_ice_penalty(q_ice):
    """约束冬季(12-2 月)冰源出流接近 0（避免非物理冬季融冰/冰源径流）。"""
    if not YEARBOOK_MODE:
        return 0.0, None
    if OBJ_W_WINTER_ICE <= 0:
        return 0.0, None
    if WINTER_MASK is None:
        return 0.0, None
    qi = np.asarray(q_ice, dtype=np.float64)
    mean_w = float(np.nanmean(qi[WINTER_MASK])) if np.any(WINTER_MASK) else float("nan")
    if (not np.isfinite(mean_w)) or mean_w <= WINTER_ICE_MAX_MEAN_M3S:
        return 0.0, mean_w
    diff = mean_w - WINTER_ICE_MAX_MEAN_M3S
    return float(OBJ_W_WINTER_ICE * diff * diff), mean_w


def _gm_volume_penalty(q_ice):
    """利用 GM_*.tif 的独立冰川融水序列，约束模拟冰源体积的量级（软惩罚）。

    说明：
    - 只在 GM 数据可用且 OBJ_W_GM_VOL>0 时启用；否则返回 0。
    - 仅对融化季（MELT_MASK: 4-9 月）做体积比约束，避免冬季接近 0 时的数值不稳定。
    - 采用 ratio 的 log-distance（相对误差对称）并在 [GM_VOL_RATIO_MIN, GM_VOL_RATIO_MAX] 内不惩罚。
    """
    if OBJ_W_GM_VOL <= 0:
        return 0.0, None
    if GM_SERIES_SIM is None:
        return 0.0, None
    if MELT_MASK is None:
        return 0.0, None

    qi = np.asarray(q_ice, dtype=np.float64)
    gm = np.asarray(GM_SERIES_SIM, dtype=np.float64)
    n = min(len(qi), len(gm), len(MELT_MASK))
    if n <= 0:
        return 0.0, None

    mask = MELT_MASK[:n]
    sim_sum = float(np.nansum(qi[:n][mask]))
    gm_sum = float(np.nansum(gm[:n][mask]))
    if (not np.isfinite(sim_sum)) or (not np.isfinite(gm_sum)) or (gm_sum <= 0.0):
        return 0.0, None

    ratio = sim_sum / (gm_sum + EPS)
    if GM_VOL_RATIO_MIN <= ratio <= GM_VOL_RATIO_MAX:
        return 0.0, float(ratio)

    if ratio > GM_VOL_RATIO_MAX:
        excess = float(np.log(ratio / GM_VOL_RATIO_MAX))
    else:
        excess = float(np.log(GM_VOL_RATIO_MIN / max(ratio, EPS)))

    return float(OBJ_W_GM_VOL * excess * excess), float(ratio)


def _compute_objective_terms(q_sim, q_rain, q_snow, q_ice):
    """Compute objective + key metrics from already-simulated series (avoid duplicate simulation)."""
    q_sim = np.asarray(q_sim, dtype=np.float64)
    q_rain = np.asarray(q_rain, dtype=np.float64)
    q_snow = np.asarray(q_snow, dtype=np.float64)
    q_ice = np.asarray(q_ice, dtype=np.float64)

    # Source closure (routed); should hold exactly when routing is linear and parameters are valid.
    resid = q_sim - (q_rain + q_snow + q_ice)
    max_abs_resid = float(np.nanmax(np.abs(resid))) if resid.size else float("nan")

    q_sim_calib = q_sim[CALIB_MASK]
    q_sim_valid = q_sim[VALID_MASK]

    min_cal = min(len(q_sim_calib), len(Q_OBS_CALIB))
    min_val = min(len(q_sim_valid), len(Q_OBS_VALID))
    obs_cal = Q_OBS_CALIB[:min_cal]
    sim_cal = q_sim_calib[:min_cal]
    obs_val = Q_OBS_VALID[:min_val]
    sim_val = q_sim_valid[:min_val]

    nse_cal = nse_numba(obs_cal, sim_cal)
    nse_val = nse_numba(obs_val, sim_val)
    nse_log_cal = nse_numba(np.log1p(np.clip(obs_cal, 0.0, None)), np.log1p(np.clip(sim_cal, 0.0, None)))
    nse_log_val = nse_numba(np.log1p(np.clip(obs_val, 0.0, None)), np.log1p(np.clip(sim_val, 0.0, None)))

    den_cal = np.nansum(obs_cal)
    pbias_cal = 0.0 if den_cal <= 0 else 100.0 * np.nansum(sim_cal - obs_cal) / den_cal

    den_val = np.nansum(obs_val)
    pbias_val = 0.0 if den_val <= 0 else 100.0 * np.nansum(sim_val - obs_val) / den_val

    obj_fit = (1.0 - nse_cal) + OBJ_W_VAL * (1.0 - nse_val)
    obj_log = (1.0 - nse_log_cal) + OBJ_W_VAL * (1.0 - nse_log_val)
    obj_base = (
        obj_fit
        + OBJ_W_LOG * obj_log
        + OBJ_W_PBIAS_CAL * (abs(pbias_cal) / 100.0)
        + OBJ_W_PBIAS_VAL * (abs(pbias_val) / 100.0)
    )

    pen_season, frac_mayoct = _seasonality_penalty(q_sim)
    pen_off, worst_excess_zero = _offseason_monthly_penalty(q_sim)
    pen_wice, winter_ice_mean = _winter_ice_penalty(q_ice)
    pen_gm, gm_ratio = _gm_volume_penalty(q_ice)
    obj = obj_base
    if np.isfinite(pen_season) and pen_season > 0.0:
        obj += pen_season
    if np.isfinite(pen_off) and pen_off > 0.0:
        obj += pen_off
    if np.isfinite(pen_wice) and pen_wice > 0.0:
        obj += pen_wice
    if np.isfinite(pen_gm) and pen_gm > 0.0:
        obj += pen_gm

    return {
        'obj': float(obj),
        'obj_base': float(obj_base),
        'obj_fit': float(obj_fit),
        'obj_log': float(obj_log),
        'pen_season': float(pen_season),
        'pen_off': float(pen_off),
        'pen_wice': float(pen_wice),
        'pen_gm': float(pen_gm),
        'mayoct_frac': frac_mayoct,
        'offseason_worst_excess_zero_frac': worst_excess_zero,
        'winter_ice_mean_m3s': winter_ice_mean,
        'gm_ratio_ice_over_gm': gm_ratio,
        'nse_cal': float(nse_cal),
        'nse_val': float(nse_val),
        'nse_log_cal': float(nse_log_cal),
        'nse_log_val': float(nse_log_val),
        'pbias_cal': float(pbias_cal),
        'pbias_val': float(pbias_val),
        'source_closure_max_abs_resid_m3s': float(max_abs_resid),
    }


def objective(x):
    """Objective function."""
    global eval_count, best_score, best_params

    try:
        eval_count += 1

        # 硬物理约束（先检查，避免无谓模拟）
        try:
            # x 顺序见 param_names：
            # 0 TT,1 FC,2 BETA,3 LP,4 RFCF,5 SFCF,6 CFR,7 CWH,8 CFMAX_mid,9 CFMAX_high,
            # 10 K,11 K1,12 K2,13 PERC,14 K_MUSK,15 ICE_FACTOR
            k0 = float(x[10])
            k1 = float(x[11])
            k2 = float(x[12])
            k_musk = float(x[14])
            if not (k0 > k1 > k2 > 0.0):
                return BAD_OBJ
            # Muskingum 系数非负域：保证分量路由线性与闭合
            if (k_musk < K_MUSK_MIN_SAFE) or (k_musk > K_MUSK_MAX_SAFE):
                return BAD_OBJ
        except Exception:
            # 若解包失败，返回劣解
            return BAD_OBJ

        q_sim, q_rain, q_snow, q_ice = run_simulation(x)
        terms = _compute_objective_terms(q_sim, q_rain, q_snow, q_ice)

        # 分量闭合检查（避免路由系数异常/截断导致的非线性破坏）
        if args.glacier_mode != "series":
            if (not np.isfinite(terms['source_closure_max_abs_resid_m3s'])) or (terms['source_closure_max_abs_resid_m3s'] > 1e-3):
                return BAD_OBJ

        if (not np.isfinite(terms['nse_cal'])) or (not np.isfinite(terms['nse_val'])) or (not np.isfinite(terms['nse_log_cal'])) or (not np.isfinite(terms['nse_log_val'])):
            return BAD_OBJ

        return terms['obj']

    except Exception as e:
        print(f"Error: {e}")
        return BAD_OBJ


def objective_with_logging(x):
    """Objective function with logging (for single-worker mode)."""
    global eval_count, best_score, best_params

    try:
        eval_count += 1

        # 硬物理约束（先检查，避免无谓模拟）
        try:
            k0 = float(x[10])
            k1 = float(x[11])
            k2 = float(x[12])
            k_musk = float(x[14])
            if not (k0 > k1 > k2 > 0.0):
                return BAD_OBJ
            if (k_musk < K_MUSK_MIN_SAFE) or (k_musk > K_MUSK_MAX_SAFE):
                return BAD_OBJ
        except Exception:
            return BAD_OBJ

        q_sim, q_rain, q_snow, q_ice = run_simulation(x)
        terms = _compute_objective_terms(q_sim, q_rain, q_snow, q_ice)

        if args.glacier_mode != "series":
            if (not np.isfinite(terms['source_closure_max_abs_resid_m3s'])) or (terms['source_closure_max_abs_resid_m3s'] > 1e-3):
                return BAD_OBJ

        nse_cal = terms['nse_cal']
        nse_val = terms['nse_val']
        nse_log_cal = terms['nse_log_cal']
        nse_log_val = terms['nse_log_val']
        pbias_cal = terms['pbias_cal']
        pbias_val = terms['pbias_val']
        obj = terms['obj']
        frac_mayoct = terms['mayoct_frac']
        pen_off = terms['pen_off']
        pen_wice = terms['pen_wice']
        pen_gm = terms.get('pen_gm', 0.0)
        wice_mean = terms['winter_ice_mean_m3s']
        gm_ratio = terms.get('gm_ratio_ice_over_gm', None)

        if (not np.isfinite(nse_cal)) or (not np.isfinite(nse_val)) or (not np.isfinite(nse_log_cal)) or (not np.isfinite(nse_log_val)):
            return BAD_OBJ

        improved = False

        # Track best by calibration NSE (主指标仍以率定期 NSE 为主)
        if nse_cal > best_score:
            best_score = nse_cal
            best_params = np.array(x).copy()
            improved = True

        if eval_count % args.log_every == 0 or improved:
            elapsed = time.time() - start_time
            rate = eval_count / elapsed if elapsed > 0 else 0
            season_txt = ""
            if frac_mayoct is not None and np.isfinite(frac_mayoct):
                season_txt = f" MayOctFrac={frac_mayoct:.3f}"
            off_txt = f" OffPen={pen_off:.3f}" if np.isfinite(pen_off) and pen_off > 0 else ""
            wice_txt = ""
            if np.isfinite(pen_wice) and pen_wice > 0 and wice_mean is not None and np.isfinite(wice_mean):
                wice_txt = f" WIceMean={wice_mean:.2f}m3s"
            gm_txt = ""
            if gm_ratio is not None and np.isfinite(gm_ratio):
                # 仅在启用 GM 约束时输出（ratio=模拟冰源体积/GM体积，融化季 4-9 月）
                gm_txt = f" GMratio={gm_ratio:.2f}"
            log_msg(
                f"[{eval_count}] NSE_cal={nse_cal:.4f} NSElog_cal={nse_log_cal:.4f} "
                f"NSE_val={nse_val:.4f} NSElog_val={nse_log_val:.4f} "
                f"PBIAS_cal={pbias_cal:+.2f}% PBIAS_val={pbias_val:+.2f}% "
                f"OBJ={obj:.4f}{season_txt}{off_txt}{wice_txt}{gm_txt} Best(NSE_cal)={best_score:.4f} Rate={rate:.1f}/s"
            )

        return obj

    except Exception as e:
        log_msg(f"Error: {e}")
        return BAD_OBJ


def make_refine_bounds(x_best, bounds, shrink=0.15, min_width=0.05):
    """Create a narrower bounds box around x_best for a second-stage DE refinement.

    说明：本项目单次模拟代价较高，SciPy 默认 polish(L-BFGS-B) 会触发大量串行有限差分评估，
    容易出现“step 80 后长时间无输出”的假卡死。用小范围二次 DE 代替 polish，既能并行，也更稳定可控。
    """
    x_best = np.asarray(x_best, dtype=np.float64)
    out = []
    for xi, (lo, hi) in zip(x_best, bounds):
        lo = float(lo)
        hi = float(hi)
        r = hi - lo
        if r <= 0:
            out.append((lo, hi))
            continue
        w = max(shrink * r, min_width * r)
        lo2 = max(lo, xi - w)
        hi2 = min(hi, xi + w)
        # 保证区间非空
        if (hi2 - lo2) < (1e-12 * max(1.0, abs(xi))):
            lo2, hi2 = lo, hi
        out.append((lo2, hi2))
    return out


def de_callback(xk, convergence):
    """每一代结束后的回调：记录进度与关键指标（在主线程执行，适用于多线程/单线程）。"""
    global gen_count, best_score, best_params

    gen_count += 1
    try:
        q_sim, q_rain, q_snow, q_ice = run_simulation(xk)
        terms = _compute_objective_terms(q_sim, q_rain, q_snow, q_ice)

        nse_cal = terms['nse_cal']
        nse_val = terms['nse_val']
        nse_log_cal = terms['nse_log_cal']
        nse_log_val = terms['nse_log_val']
        pbias_cal = terms['pbias_cal']
        pbias_val = terms['pbias_val']
        obj = terms['obj']
        frac_mayoct = terms['mayoct_frac']
        pen_off = terms['pen_off']
        pen_wice = terms['pen_wice']
        wice_mean = terms['winter_ice_mean_m3s']
        gm_ratio = terms.get('gm_ratio_ice_over_gm', None)

        if args.glacier_mode != "series":
            if (not np.isfinite(terms['source_closure_max_abs_resid_m3s'])) or (terms['source_closure_max_abs_resid_m3s'] > 1e-3):
                obj = BAD_OBJ

        if np.isfinite(nse_cal) and (nse_cal > best_score):
            best_score = nse_cal
            best_params = np.array(xk).copy()

        append_progress(
            gen_count,
            float(nse_cal), float(nse_log_cal),
            float(nse_val), float(nse_log_val),
            float(pbias_cal), float(pbias_val),
            float(obj), float(convergence),
        )

        stage = "REFINE" if (PROGRESS_FILE_REFINE and PROGRESS_FILE == PROGRESS_FILE_REFINE) else "GLOBAL"
        season_txt = ""
        if frac_mayoct is not None and np.isfinite(frac_mayoct):
            season_txt = f" MayOctFrac={frac_mayoct:.3f}"
        off_txt = f" OffPen={pen_off:.3f}" if np.isfinite(pen_off) and pen_off > 0 else ""
        wice_txt = ""
        if np.isfinite(pen_wice) and pen_wice > 0 and wice_mean is not None and np.isfinite(wice_mean):
            wice_txt = f" WIceMean={wice_mean:.2f}m3s"
        gm_txt = ""
        if gm_ratio is not None and np.isfinite(gm_ratio):
            gm_txt = f" GMratio={gm_ratio:.2f}"
        log_msg(
            f"[{stage} Gen {gen_count}] NSE_cal={nse_cal:.4f} NSElog_cal={nse_log_cal:.4f} "
            f"NSE_val={nse_val:.4f} NSElog_val={nse_log_val:.4f} "
            f"PBIAS_cal={pbias_cal:+.2f}% PBIAS_val={pbias_val:+.2f}% "
            f"OBJ={obj:.4f}{season_txt}{off_txt}{wice_txt}{gm_txt} conv={convergence:.3e}"
        )

    except Exception as e:
        log_msg(f"[Gen {gen_count}] Callback error: {e}")

    # 返回 True 可提前停止优化；这里始终继续
    return False


def _calc_glacier_area_km2():
    """Compute glacier area used by the model grid (km2, already area-corrected to catchment area)."""
    if GLACIER_MASK is None or PX_AREA is None or FLOW_ACC is None:
        return 0.0
    valid_mask = ~np.isnan(FLOW_ACC)
    return float(np.sum(PX_AREA[valid_mask & GLACIER_MASK]) * AREA_COEF)


def _corr(a, b, mask=None):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if mask is None:
        m = np.isfinite(a) & np.isfinite(b)
    else:
        mask = np.asarray(mask, dtype=bool)
        m = mask & np.isfinite(a) & np.isfinite(b)
    if m.sum() < 2:
        return float("nan")
    return float(np.corrcoef(a[m], b[m])[0, 1])


def _center_of_mass_doy(dates, x):
    x = np.asarray(x, dtype=np.float64)
    if np.nansum(x) <= 0:
        return float("nan")
    doy = pd.to_datetime(dates).dayofyear.values.astype(np.float64)
    return float(np.nansum(doy * x) / np.nansum(x))


def _sanitize_for_json(obj):
    """Make dict/list/json-scalars safe (replace NaN/Inf with None)."""
    if isinstance(obj, dict):
        return {k: _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_for_json(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    return obj


def compute_and_save_diagnostics(run_dir, dates, q_sim, q_rain, q_snow, q_ice,
                                 q_snow_glacier, q_glacier_total,
                                 q_ice_raw, q_snow_glacier_raw, skip_gm=False):
    """Compute physical self-checks + independent GM comparison and save to run directory."""
    diag = {}

    # 1) Source closure (routed)
    resid = q_sim - (q_rain + q_snow + q_ice)
    diag['source_closure'] = {
        'max_abs_residual_m3s': float(np.nanmax(np.abs(resid))),
        'mean_abs_residual_m3s': float(np.nanmean(np.abs(resid))),
    }

    # 1.5) Component fractions (routed, full period; glacier_total = q_ice + q_snow_glacier)
    try:
        snow_nogl = np.clip(np.asarray(q_snow, dtype=np.float64) - np.asarray(q_snow_glacier, dtype=np.float64), 0.0, None)
        glacier_total = np.asarray(q_glacier_total, dtype=np.float64)
        den_q = float(np.nansum(np.asarray(q_sim, dtype=np.float64)))
        if den_q > 0:
            diag['component_fractions_routed'] = {
                'rain': float(np.nansum(np.asarray(q_rain, dtype=np.float64)) / den_q),
                'snow_noglacier': float(np.nansum(snow_nogl) / den_q),
                'glacier_total': float(np.nansum(glacier_total) / den_q),
                'ice': float(np.nansum(np.asarray(q_ice, dtype=np.float64)) / den_q),
                'snow_glacier': float(np.nansum(np.asarray(q_snow_glacier, dtype=np.float64)) / den_q),
            }
        else:
            diag['component_fractions_routed'] = None
    except Exception:
        diag['component_fractions_routed'] = None

    # 2) Glacier area + ice intensity checks (use raw series, closer to GM definition)
    glacier_area_km2 = _calc_glacier_area_km2()
    diag['glacier_area_km2'] = round(glacier_area_km2, 3)

    months = pd.to_datetime(dates).month.values.astype(int)
    melt_mask = (months >= 4) & (months <= 9)
    winter_mask = (months == 12) | (months <= 2)

    if glacier_area_km2 > 0:
        ice_mm_d = (np.asarray(q_ice_raw, dtype=np.float64) * 86.4) / glacier_area_km2
        diag['ice_intensity_over_glacier_mm_d'] = {
            'max': float(np.nanmax(ice_mm_d)),
            'p99': float(np.nanpercentile(ice_mm_d, 99)),
            'mean_aprsep': float(np.nanmean(ice_mm_d[melt_mask])),
        }

        df_tmp = pd.DataFrame({'date': pd.to_datetime(dates), 'q_ice_raw': q_ice_raw})
        annual = df_tmp.groupby(df_tmp['date'].dt.year)['q_ice_raw'].sum().astype(np.float64)
        annual_mm = (annual.values * 86.4) / glacier_area_km2
        diag['ice_depth_over_glacier_mm_yr'] = {
            'min': float(np.nanmin(annual_mm)) if annual_mm.size else float("nan"),
            'p50': float(np.nanmedian(annual_mm)) if annual_mm.size else float("nan"),
            'p90': float(np.nanpercentile(annual_mm, 90)) if annual_mm.size else float("nan"),
            'max': float(np.nanmax(annual_mm)) if annual_mm.size else float("nan"),
        }
    else:
        diag['ice_intensity_over_glacier_mm_d'] = None
        diag['ice_depth_over_glacier_mm_yr'] = None

    diag['winter_ice'] = {
        'mean_m3s': float(np.nanmean(np.asarray(q_ice, dtype=np.float64)[winter_mask])),
        'max_m3s': float(np.nanmax(np.asarray(q_ice, dtype=np.float64)[winter_mask])),
    }
    q_sim_arr = np.asarray(q_sim, dtype=np.float64)
    if np.any(winter_mask):
        wq = q_sim_arr[winter_mask]
        diag['winter_qsim'] = {
            'mean_m3s': float(np.nanmean(wq)),
            'p95_m3s': float(np.nanpercentile(wq, 95)),
            'max_m3s': float(np.nanmax(wq)),
        }
    else:
        diag['winter_qsim'] = None

    # 3) Seasonality summary
    df_m = pd.DataFrame({
        'date': pd.to_datetime(dates),
        'q_sim': q_sim,
        'q_rain': q_rain,
        'q_snow': q_snow,
        'q_ice': q_ice,
        'q_snow_glacier': q_snow_glacier,
        'q_glacier_total': q_glacier_total,
    })
    monthly_mean = df_m.groupby(df_m['date'].dt.month)[['q_sim', 'q_rain', 'q_snow', 'q_ice', 'q_snow_glacier', 'q_glacier_total']].mean()
    diag['peak_month'] = {k: int(monthly_mean[k].idxmax()) for k in monthly_mean.columns}
    diag['center_of_mass_doy'] = {
        'q_snow': round(_center_of_mass_doy(dates, q_snow), 2),
        'q_ice': round(_center_of_mass_doy(dates, q_ice), 2),
        'q_glacier_total': round(_center_of_mass_doy(dates, q_glacier_total), 2),
    }
    # 补充：全期/逐年季节性占比（用于检查“冬季是否异常高流量”）
    try:
        qsim_tot = float(np.nansum(q_sim_arr))
        if qsim_tot > 0 and SIM_MONTH is not None and len(SIM_MONTH) == len(q_sim_arr):
            mayoct_mask2 = (SIM_MONTH >= 5) & (SIM_MONTH <= 10)
            novapr_mask2 = (SIM_MONTH >= 11) | (SIM_MONTH <= 4)
            djf_mask2 = (SIM_MONTH == 12) | (SIM_MONTH <= 2)
            frac_mayoct_total = float(np.nansum(q_sim_arr[mayoct_mask2]) / (qsim_tot + EPS))
            frac_novapr_total = float(np.nansum(q_sim_arr[novapr_mask2]) / (qsim_tot + EPS))
            frac_djf_total = float(np.nansum(q_sim_arr[djf_mask2]) / (qsim_tot + EPS))

            # 逐年 5-10 月占比
            df_frac = pd.DataFrame({'date': pd.to_datetime(dates), 'q': q_sim_arr})
            df_frac['year'] = df_frac['date'].dt.year
            df_frac['month'] = df_frac['date'].dt.month
            frac_list = []
            for _, g in df_frac.groupby('year'):
                tot_y = float(g['q'].sum())
                if tot_y <= 0:
                    continue
                mayoct_y = float(g.loc[g['month'].between(5, 10), 'q'].sum())
                frac_list.append(mayoct_y / tot_y)
            frac_arr = np.asarray(frac_list, dtype=np.float64)
            diag['seasonality_qsim'] = {
                'frac_mayoct_total': frac_mayoct_total,
                'frac_novapr_total': frac_novapr_total,
                'frac_djf_total': frac_djf_total,
                'frac_mayoct_year_mean': float(np.nanmean(frac_arr)) if frac_arr.size else None,
                'frac_mayoct_year_std': float(np.nanstd(frac_arr)) if frac_arr.size else None,
                'frac_mayoct_year_min': float(np.nanmin(frac_arr)) if frac_arr.size else None,
                'frac_mayoct_year_max': float(np.nanmax(frac_arr)) if frac_arr.size else None,
                'yearbook_target_frac_mayoct': (
                    float(YEARBOOK_MAYOCT_FRAC)
                    if (YEARBOOK_MAYOCT_FRAC is not None and np.isfinite(YEARBOOK_MAYOCT_FRAC))
                    else (YEARBOOK_MAYOCT_FRAC_DEFAULT if YEARBOOK_MODE else None)
                ),
            }
        else:
            diag['seasonality_qsim'] = None
    except Exception:
        diag['seasonality_qsim'] = None

    # 4) Independent GM comparison (glacier melt rasters)
    gm = None
    if (not skip_gm) and os.path.isdir(GLACIER_MELT_DIR):
        try:
            print("  [DIAG] Loading GM_*.tif for independent glacier-melt comparison (one-time)...")
            gm = load_glacier_melt_series(pd.to_datetime(dates).min().date(), pd.to_datetime(dates).max().date())
        except Exception as e:
            gm = None
            diag['gm_error'] = str(e)

    if gm is not None and len(gm) == len(dates):
        q_glacier_total_raw = np.asarray(q_snow_glacier_raw, dtype=np.float64) + np.asarray(q_ice_raw, dtype=np.float64)
        diag['gm_comparison'] = {
            'corr_aprsep_gm_vs_ice': _corr(gm, q_ice_raw, melt_mask),
            'corr_aprsep_gm_vs_glacier_total': _corr(gm, q_glacier_total_raw, melt_mask),
            'sum_ratio_ice_over_gm': float(np.nansum(q_ice_raw) / np.nansum(gm)) if np.nansum(gm) > 0 else float("nan"),
            'sum_ratio_glacier_total_over_gm': float(np.nansum(q_glacier_total_raw) / np.nansum(gm)) if np.nansum(gm) > 0 else float("nan"),
            'gm_mean_m3s': float(np.nanmean(gm)),
            'gm_max_m3s': float(np.nanmax(gm)),
        }
        # Save a compact time series for cryosphere diagnosis
        df_cryo = pd.DataFrame({
            'date': pd.to_datetime(dates),
            'gm_melt_m3s': gm,
            'q_ice_raw_m3s': q_ice_raw,
            'q_snow_glacier_raw_m3s': q_snow_glacier_raw,
            'q_glacier_total_raw_m3s': q_glacier_total_raw,
        })
        df_cryo.to_csv(os.path.join(run_dir, 'cryo_diagnostics.csv'), index=False)
    else:
        diag['gm_comparison'] = None

    # 5) Hard-check flags (conservative thresholds; only for warnings, not used to optimize)
    checks = {}
    checks['source_closure_ok'] = diag['source_closure']['max_abs_residual_m3s'] < 1e-4
    if glacier_area_km2 > 0 and diag['ice_intensity_over_glacier_mm_d'] is not None:
        checks['ice_mm_d_p99_ok'] = diag['ice_intensity_over_glacier_mm_d']['p99'] < 60.0
        checks['ice_mm_yr_max_ok'] = diag['ice_depth_over_glacier_mm_yr']['max'] < 3000.0
    else:
        checks['ice_mm_d_p99_ok'] = True
        checks['ice_mm_yr_max_ok'] = True
    checks['winter_ice_small_ok'] = diag['winter_ice']['mean_m3s'] < WINTER_ICE_MAX_MEAN_M3S
    try:
        if YEARBOOK_MODE and diag.get('seasonality_qsim') and diag['seasonality_qsim'].get('frac_mayoct_total') is not None:
            target = diag['seasonality_qsim'].get('yearbook_target_frac_mayoct', YEARBOOK_MAYOCT_FRAC_DEFAULT)
            min_frac = max(0.0, float(target) - YEARBOOK_MAYOCT_FRAC_SLACK)
            checks['seasonality_mayoct_ok'] = float(diag['seasonality_qsim']['frac_mayoct_total']) >= min_frac
        else:
            checks['seasonality_mayoct_ok'] = True
    except Exception:
        checks['seasonality_mayoct_ok'] = True
    diag['checks'] = checks

    # Save diagnostics JSON
    diag_file = os.path.join(run_dir, 'diagnostics.json')
    diag_clean = _sanitize_for_json(diag)
    with open(diag_file, 'w', encoding='utf-8') as f:
        json.dump(diag_clean, f, indent=2, ensure_ascii=False, allow_nan=False)

    return diag_clean


def main():
    global args

    args = parse_args()
    # Default: always do GM comparison when available (only diagnostics, not used in optimization).
    # Can be disabled via --skip-gm.
    if args.skip_gm is None:
        args.skip_gm = False
    setup_logging()

    global PREC_DIR
    if args.prec_dir:
        PREC_DIR = args.prec_dir
    else:
        PREC_DIR = PREC_DIR_CMFD if args.prec_source == "cmfd" else PREC_DIR_MSWEP

    # Analysis mode: read an existing run directory and output diagnostics without recalibration
    if args.analyze_run:
        run_dir = args.analyze_run
        meta_path = os.path.join(run_dir, 'metadata.json')
        if not os.path.exists(meta_path):
            raise FileNotFoundError(f"--analyze-run provided but metadata.json not found: {meta_path}")
        with open(meta_path, 'r', encoding='utf-8') as f:
            meta = json.load(f)

        # Override runtime config from metadata for reproducibility
        args.prec_source = meta.get('data_sources', {}).get('prec_source', args.prec_source)
        args.glacier_mode = meta.get('data_sources', {}).get('glacier_mode', args.glacier_mode)
        if (not args.prec_dir) and meta.get('data_sources', {}).get('prec_dir'):
            PREC_DIR = meta['data_sources']['prec_dir']
        else:
            PREC_DIR = PREC_DIR_CMFD if args.prec_source == "cmfd" else PREC_DIR_MSWEP

        # Load optimized parameters (ordered)
        if 'optimized_params' not in meta:
            raise KeyError(f"metadata.json missing optimized_params: {meta_path}")
        x = [float(meta['optimized_params'][name]) for name in param_names]

        print("="*70)
        print("="*70)
        print(f"[ANALYZE] 目标流域 HBV V6 运行结果诊断: {run_dir}")
        print("="*70)
        print("="*70)
        print(f"Source: {args.prec_source} | Glacier mode: {args.glacier_mode}")
        print(f"Using parameters from: {meta_path}")

        print("\n" + "-"*70)
        warmup_jit()

        load_all_data()

        (q_sim, q_rain, q_snow, q_ice,
         q_snow_glacier, q_glacier_total,
         q_total_raw, q_rain_raw, q_snow_raw, q_ice_raw, q_snow_glacier_raw) = run_simulation_detailed(x)

        # Save detailed simulation (keep original simulation.csv untouched)
        min_len = min(len(q_sim), len(Q_OBS_FULL))
        df_det = pd.DataFrame({
            'date': SIM_DATES[:min_len],
            'q_sim': q_sim[:min_len],
            'q_obs': Q_OBS_FULL[:min_len],
            'q_obs_obj': (Q_OBS_OBJ_FULL[:min_len] if Q_OBS_OBJ_FULL is not None else np.full(min_len, np.nan)),
            'obs_masked': (OBS_OBJ_MASK[:min_len].astype(np.int8) if OBS_OBJ_MASK is not None else np.zeros(min_len, dtype=np.int8)),
            'q_rain': q_rain[:min_len],
            'q_snow': q_snow[:min_len],
            'q_snow_noglacier': np.clip(q_snow[:min_len] - q_snow_glacier[:min_len], 0.0, None),
            'q_ice': q_ice[:min_len],
            'q_snow_glacier': q_snow_glacier[:min_len],
            'q_glacier_total': q_glacier_total[:min_len],
            'q_total_raw': q_total_raw[:min_len],
            'q_rain_raw': q_rain_raw[:min_len],
            'q_snow_raw': q_snow_raw[:min_len],
            'q_snow_noglacier_raw': np.clip(q_snow_raw[:min_len] - q_snow_glacier_raw[:min_len], 0.0, None),
            'q_ice_raw': q_ice_raw[:min_len],
            'q_snow_glacier_raw': q_snow_glacier_raw[:min_len],
            'q_glacier_total_raw': (q_ice_raw[:min_len] + q_snow_glacier_raw[:min_len]),
        })
        det_csv = os.path.join(run_dir, 'simulation_detailed.csv')
        df_det.to_csv(det_csv, index=False)

        diag = compute_and_save_diagnostics(
            run_dir,
            SIM_DATES[:min_len],
            q_sim[:min_len], q_rain[:min_len], q_snow[:min_len], q_ice[:min_len],
            q_snow_glacier[:min_len], q_glacier_total[:min_len],
            q_ice_raw[:min_len], q_snow_glacier_raw[:min_len],
            skip_gm=args.skip_gm,
        )

        print("\n" + "=" * 60)
        print("DIAGNOSTICS SAVED")
        print("=" * 60)
        print(f"  - simulation_detailed.csv: {det_csv}")
        print(f"  - diagnostics.json       : {os.path.join(run_dir, 'diagnostics.json')}")
        if os.path.exists(os.path.join(run_dir, 'cryo_diagnostics.csv')):
            print(f"  - cryo_diagnostics.csv   : {os.path.join(run_dir, 'cryo_diagnostics.csv')}")
        print("Checks:", diag.get('checks', {}))
        return

    print("="*70)
    print("="*70)
    print(f"目标流域 HBV 半分布式模型 - V6 空间分区率定 ({len(param_names)}参数)")
    print("="*70)
    print("="*70)
    print(f"\nRun ID: {RUN_ID}")
    print(f"Workers: {args.workers}")
    print(f"Max iterations: {args.maxiter}")
    print(f"Population: {args.popsize * len(PARAM_BOUNDS)}")
    print(f"Precip source: {args.prec_source}")
    print(f"Precip dir: {PREC_DIR}")
    print(f"Glacier mode: {args.glacier_mode}")
    print(f"\n优化参数: {param_names}")
    print("模型设定:")
    print(f"  - 冰川模式: {args.glacier_mode} (SWE阈值={SWE_ICE_THRESHOLD_MM} mm w.e.)")
    print("  - PET: ERA5 FAO-56 Penman-Monteith, ALPHA/ECORR 固定为 0")
    print("  - CFMAX: 中/高海拔分区率定")
    print("  - 物理约束: K > K1 > K2; Muskingum K_MUSK 位于系数非负安全域")
    print("  - 时间划分: 06-08预热, 09-17率定, 18-20验证")

    # Warmup JIT
    print("\n" + "-"*70)
    warmup_jit()

    # Load data
    if args.quick_test:
        print(f"\n[INFO] Quick test mode: loading only {args.quick_days} days...")
        quick_end = (datetime.strptime(WARMUP_START, "%Y-%m-%d") + timedelta(days=max(args.quick_days, 1) - 1)).strftime("%Y-%m-%d")
        load_all_data(end_date_override=quick_end, skip_obs=True)
    else:
        load_all_data()

    # Benchmark
    print("\n" + "-"*70)
    print("Performance test:")
    # V6 benchmark params (must match `param_names` order):
    # TT, FC, BETA, LP, RFCF, SFCF, CFR, CWH, CFMAX_mid, CFMAX_high, K, K1, K2, PERC, K_MUSK, ICE_FACTOR
    x_test = [-1.2, 1200.0, 1.0, 0.95, 1.0, 1.05, 0.05, 0.05, 3.5, 5.5, 0.25, 0.05, 0.005, 1.8, 2.5, 2.0]
    t0 = time.time()
    q_test, _, _, _ = run_simulation(x_test)
    t_single = time.time() - t0
    print(f"  Single evaluation: {t_single:.3f} seconds")

    if args.quick_test:
        print("\n[INFO] Quick test completed successfully.")
        return

    n_evals = args.popsize * len(PARAM_BOUNDS) * (args.maxiter + 1)
    est_hours_serial = n_evals * t_single / 3600
    print(f"  Estimated evaluations: ~{n_evals}")
    print(f"  Estimated time (serial): ~{est_hours_serial:.2f} hours")
    if args.workers and args.workers > 1:
        est_hours_wall = est_hours_serial / max(int(args.workers), 1)
        print(f"  Estimated time (wall, ~{args.workers} workers): ~{est_hours_wall:.2f} hours (rough)")

    # Run optimization
    print("\n" + "="*70)
    print("Starting optimization...")
    print("="*70)
    if args.polish:
        print("[INFO] polish=True: 完成 maxiter 代进化后将执行局部寻优(polish/L-BFGS-B)。该阶段可能耗时较长且无 step 输出。")
    else:
        print("[INFO] 默认 polish=False：避免高代价串行有限差分导致的“假卡死”。改用二次 DE 小范围精修（并行、可输出进度）。")

    # Stage 1: global DE
    result_global = None
    result_refine = None
    # Stage 2 defaults (only used when polish=False and workers>1)
    # Refine stage scales with global maxiter (cap to keep runtime predictable)
    REFINE_MAXITER = min(20, max(6, int(round(args.maxiter * 0.25))))
    REFINE_POPSIZE = max(6, int(round(args.popsize * 0.5)))
    # refine 阶段应更严格（更小 tol），否则容易 1-2 代提前停止
    REFINE_TOL = max(args.tol * 0.5, 1e-4)
    REFINE_SHRINK = 0.15
    REFINE_MIN_WIDTH = 0.05
    do_refine = (not args.polish) and (args.workers > 1)

    if args.workers > 1:
        # Multi-worker mode: use ThreadPool (Win compatible, no pickle overhead)
        print(f"  Using {args.workers} threads (Win-safe)...")

        # We can pass the map function of a ThreadPool directly to 'workers'
        # Note: ThreadPool shares memory, so no need for initializer/SHARED_DATA
        with ThreadPool(processes=args.workers) as pool:
            result_global = differential_evolution(
                objective,
                PARAM_BOUNDS,
                strategy='best1bin',
                maxiter=args.maxiter,
                popsize=args.popsize,
                seed=args.seed,
                tol=args.tol,
                disp=True,
                polish=args.polish,
                callback=de_callback,
                workers=pool.map,
                updating='deferred',
            )

            # Stage 2: refinement DE (recommended when polish=False)
            if do_refine:
                global PROGRESS_FILE, PROGRESS_FILE_REFINE, PROGRESS_FILE_MAIN, gen_count
                PROGRESS_FILE_REFINE = os.path.join(LOG_DIR, f"v6_progress_refine_{RUN_ID}.csv")
                with open(PROGRESS_FILE_REFINE, 'w') as f:
                    f.write("timestamp,gen,nse_cal,nse_log_cal,nse_val,nse_log_val,pbias_cal,pbias_val,obj,convergence,elapsed_sec\n")
                PROGRESS_FILE = PROGRESS_FILE_REFINE
                gen_count = 0

                bounds_refine = make_refine_bounds(
                    result_global.x, PARAM_BOUNDS,
                    shrink=REFINE_SHRINK, min_width=REFINE_MIN_WIDTH
                )
                print("\n" + "-"*70)
                print("二次精修（Stage 2/2）：缩小边界的 DE（并行/有输出，替代 polish）")
                print(f"  refine_maxiter={REFINE_MAXITER} refine_popsize={REFINE_POPSIZE} refine_tol={REFINE_TOL}")
                print(f"  shrink={REFINE_SHRINK} min_width={REFINE_MIN_WIDTH}")
                print("-"*70)

                result_refine = differential_evolution(
                    objective,
                    bounds_refine,
                    strategy='best1bin',
                    maxiter=REFINE_MAXITER,
                    popsize=REFINE_POPSIZE,
                    seed=args.seed + 1,
                    tol=REFINE_TOL,
                    disp=True,
                    polish=False,
                    callback=de_callback,
                    workers=pool.map,
                    updating='deferred',
                )

                # Restore progress file pointer (for any later logging)
                PROGRESS_FILE = PROGRESS_FILE_MAIN
    else:
        # Single-worker mode: use logging objective
        result_global = differential_evolution(
            objective_with_logging,
            PARAM_BOUNDS,
            strategy='best1bin',
            maxiter=args.maxiter,
            popsize=args.popsize,
            seed=args.seed,
            tol=args.tol,
            disp=True,
            polish=args.polish,
            callback=de_callback,
            workers=1,
            updating='immediate',
        )
        result_refine = None

    # Pick best result by objective value (lower is better)
    result = result_global
    selected_stage = "global"
    if result_refine is not None and np.isfinite(getattr(result_refine, "fun", np.inf)):
        if getattr(result_refine, "fun", np.inf) < getattr(result_global, "fun", np.inf):
            result = result_refine
            selected_stage = "refine"

    # Results
    elapsed_total = time.time() - start_time
    nfev_global = int(getattr(result_global, "nfev", 0)) if result_global is not None else 0
    nfev_refine = int(getattr(result_refine, "nfev", 0)) if result_refine is not None else 0
    nfev = nfev_global + nfev_refine if (nfev_global + nfev_refine) > 0 else getattr(result, "nfev", eval_count)
    print("\n" + "="*70)
    print("OPTIMIZATION COMPLETE")
    print("="*70)
    print(f"Selected stage: {selected_stage}")
    print(f"\nTotal time: {elapsed_total/3600:.2f} hours ({elapsed_total:.0f} seconds)")
    print(f"Evaluations: {nfev}")
    print(f"Rate: {nfev/elapsed_total:.1f} eval/s")

    print("\nOptimal parameters:")
    for name, val in zip(param_names, result.x):
        print(f"  {name:10s} = {val:.4f}")

    # Final metrics (calibration + validation) + cryosphere diagnostics (forward run once)
    (q_sim, q_rain, q_snow, q_ice,
     q_snow_glacier, q_glacier_total,
     q_total_raw, q_rain_raw, q_snow_raw, q_ice_raw, q_snow_glacier_raw) = run_simulation_detailed(result.x)
    q_sim_calib = q_sim[CALIB_MASK]
    q_sim_valid = q_sim[VALID_MASK]

    min_calib = min(len(q_sim_calib), len(Q_OBS_CALIB))
    min_valid = min(len(q_sim_valid), len(Q_OBS_VALID))

    nse_calib = nse_numba(Q_OBS_CALIB[:min_calib], q_sim_calib[:min_calib])
    nse_valid = nse_numba(Q_OBS_VALID[:min_valid], q_sim_valid[:min_valid])
    den_cal = np.nansum(Q_OBS_CALIB[:min_calib])
    den_val = np.nansum(Q_OBS_VALID[:min_valid])
    pbias_calib = 0.0 if den_cal <= 0 else 100.0 * np.nansum(q_sim_calib[:min_calib] - Q_OBS_CALIB[:min_calib]) / den_cal
    pbias_valid = 0.0 if den_val <= 0 else 100.0 * np.nansum(q_sim_valid[:min_valid] - Q_OBS_VALID[:min_valid]) / den_val
    rmse_calib = np.sqrt(np.nanmean((Q_OBS_CALIB[:min_calib] - q_sim_calib[:min_calib])**2))
    rmse_valid = np.sqrt(np.nanmean((Q_OBS_VALID[:min_valid] - q_sim_valid[:min_valid])**2))

    # Also report metrics against raw observations (same口径 as旧版本：包含冬季0)
    obs_calib_raw = Q_OBS_FULL[CALIB_MASK]
    obs_valid_raw = Q_OBS_FULL[VALID_MASK]
    min_calib_raw = min(len(q_sim_calib), len(obs_calib_raw))
    min_valid_raw = min(len(q_sim_valid), len(obs_valid_raw))
    nse_calib_raw = nse_numba(obs_calib_raw[:min_calib_raw], q_sim_calib[:min_calib_raw])
    nse_valid_raw = nse_numba(obs_valid_raw[:min_valid_raw], q_sim_valid[:min_valid_raw])
    den_cal_raw = np.nansum(obs_calib_raw[:min_calib_raw])
    den_val_raw = np.nansum(obs_valid_raw[:min_valid_raw])
    pbias_calib_raw = 0.0 if den_cal_raw <= 0 else 100.0 * np.nansum(q_sim_calib[:min_calib_raw] - obs_calib_raw[:min_calib_raw]) / den_cal_raw
    pbias_valid_raw = 0.0 if den_val_raw <= 0 else 100.0 * np.nansum(q_sim_valid[:min_valid_raw] - obs_valid_raw[:min_valid_raw]) / den_val_raw
    rmse_calib_raw = np.sqrt(np.nanmean((obs_calib_raw[:min_calib_raw] - q_sim_calib[:min_calib_raw])**2))
    rmse_valid_raw = np.sqrt(np.nanmean((obs_valid_raw[:min_valid_raw] - q_sim_valid[:min_valid_raw])**2))

    # 组分占比（用于与文献/物理认知对照）
    # - Rain: 降雨产流
    # - Snowmelt (non-glacier): 非冰川区融雪产流（q_snow - q_snow_glacier）
    # - Glacier total: 冰川区总产流（裸冰融 + 冰川区融雪）
    rain_calib = q_rain[CALIB_MASK]
    snow_calib = q_snow[CALIB_MASK]
    ice_calib = q_ice[CALIB_MASK]
    snow_glacier_calib = q_snow_glacier[CALIB_MASK]
    snow_nogl_calib = np.clip(snow_calib - snow_glacier_calib, 0.0, None)
    glacier_calib = snow_glacier_calib + ice_calib

    rain_valid = q_rain[VALID_MASK]
    snow_valid = q_snow[VALID_MASK]
    ice_valid = q_ice[VALID_MASK]
    snow_glacier_valid = q_snow_glacier[VALID_MASK]
    snow_nogl_valid = np.clip(snow_valid - snow_glacier_valid, 0.0, None)
    glacier_valid = snow_glacier_valid + ice_valid

    den_q_cal = np.nansum(q_sim_calib)
    den_q_val = np.nansum(q_sim_valid)
    rain_frac_calib = np.nansum(rain_calib) / den_q_cal if den_q_cal > 0 else 0.0
    snow_nogl_frac_calib = np.nansum(snow_nogl_calib) / den_q_cal if den_q_cal > 0 else 0.0
    glacier_frac_calib = np.nansum(glacier_calib) / den_q_cal if den_q_cal > 0 else 0.0
    ice_frac_calib = np.nansum(ice_calib) / den_q_cal if den_q_cal > 0 else 0.0
    snow_glacier_frac_calib = np.nansum(snow_glacier_calib) / den_q_cal if den_q_cal > 0 else 0.0

    rain_frac_valid = np.nansum(rain_valid) / den_q_val if den_q_val > 0 else 0.0
    snow_nogl_frac_valid = np.nansum(snow_nogl_valid) / den_q_val if den_q_val > 0 else 0.0
    glacier_frac_valid = np.nansum(glacier_valid) / den_q_val if den_q_val > 0 else 0.0
    ice_frac_valid = np.nansum(ice_valid) / den_q_val if den_q_val > 0 else 0.0
    snow_glacier_frac_valid = np.nansum(snow_glacier_valid) / den_q_val if den_q_val > 0 else 0.0

    print(f"\nFinal metrics (calibration):")
    print(f"  NSE (objective):   {nse_calib:.4f}")
    print(f"  NSE (full):        {nse_calib_raw:.4f}")
    print(f"  PBIAS (objective): {pbias_calib:.2f}%")
    print(f"  PBIAS (full):      {pbias_calib_raw:.2f}%")
    print(f"  RMSE (objective):  {rmse_calib:.2f} m3/s")
    print(f"  RMSE (full):       {rmse_calib_raw:.2f} m3/s")
    print(f"  Rain fraction: {rain_frac_calib*100:.1f}%")
    print(f"  Snowmelt (non-glacier) fraction: {snow_nogl_frac_calib*100:.1f}%")
    print(f"  Glacier fraction: {glacier_frac_calib*100:.1f}% (ice {ice_frac_calib*100:.1f}%, glacier-snow {snow_glacier_frac_calib*100:.1f}%)")
    print(f"\nFinal metrics (validation):")
    print(f"  NSE (objective):   {nse_valid:.4f}")
    print(f"  NSE (full):        {nse_valid_raw:.4f}")
    print(f"  PBIAS (objective): {pbias_valid:.2f}%")
    print(f"  PBIAS (full):      {pbias_valid_raw:.2f}%")
    print(f"  RMSE (obj):  {rmse_valid:.2f} m3/s")
    print(f"  RMSE (raw):  {rmse_valid_raw:.2f} m3/s")
    print(f"  Rain fraction: {rain_frac_valid*100:.1f}%")
    print(f"  Snowmelt (non-glacier) fraction: {snow_nogl_frac_valid*100:.1f}%")
    print(f"  Glacier fraction: {glacier_frac_valid*100:.1f}% (ice {ice_frac_valid*100:.1f}%, glacier-snow {snow_glacier_frac_valid*100:.1f}%)")

    # Save results to dedicated run directory
    run_dir = os.path.join(PROJECT_ROOT, 'results', 'runs', f'v6_{args.prec_source}_{args.glacier_mode}_{RUN_ID}')
    os.makedirs(run_dir, exist_ok=True)

    # 1. Save simulation results CSV
    dates_sim = SIM_DATES
    min_len = min(len(q_sim), len(Q_OBS_FULL))
    df = pd.DataFrame({
        'date': dates_sim[:min_len],
        'q_sim': q_sim[:min_len],
        'q_obs': Q_OBS_FULL[:min_len],
        'q_rain': q_rain[:min_len],
        'q_snow': q_snow[:min_len],
        'q_ice': q_ice[:min_len]
    })
    result_file = os.path.join(run_dir, 'simulation.csv')
    df.to_csv(result_file, index=False)

    # 1.1 Save detailed simulation for cryosphere diagnosis (do not change simulation.csv schema)
    df_det = pd.DataFrame({
        'date': dates_sim[:min_len],
        'q_sim': q_sim[:min_len],
        'q_obs': Q_OBS_FULL[:min_len],
        'q_obs_obj': (Q_OBS_OBJ_FULL[:min_len] if Q_OBS_OBJ_FULL is not None else np.full(min_len, np.nan)),
        'obs_masked': (OBS_OBJ_MASK[:min_len].astype(np.int8) if OBS_OBJ_MASK is not None else np.zeros(min_len, dtype=np.int8)),
        'q_rain': q_rain[:min_len],
        'q_snow': q_snow[:min_len],
        'q_snow_noglacier': np.clip(q_snow[:min_len] - q_snow_glacier[:min_len], 0.0, None),
        'q_ice': q_ice[:min_len],
        'q_snow_glacier': q_snow_glacier[:min_len],
        'q_glacier_total': q_glacier_total[:min_len],
        'q_total_raw': q_total_raw[:min_len],
        'q_rain_raw': q_rain_raw[:min_len],
        'q_snow_raw': q_snow_raw[:min_len],
        'q_snow_noglacier_raw': np.clip(q_snow_raw[:min_len] - q_snow_glacier_raw[:min_len], 0.0, None),
        'q_ice_raw': q_ice_raw[:min_len],
        'q_snow_glacier_raw': q_snow_glacier_raw[:min_len],
        'q_glacier_total_raw': (q_ice_raw[:min_len] + q_snow_glacier_raw[:min_len]),
    })
    df_det.to_csv(os.path.join(run_dir, 'simulation_detailed.csv'), index=False)

    # 1.2 Physical self-checks + independent GM comparison (writes diagnostics.json / cryo_diagnostics.csv)
    diag = compute_and_save_diagnostics(
        run_dir,
        dates_sim[:min_len],
        q_sim[:min_len], q_rain[:min_len], q_snow[:min_len], q_ice[:min_len],
        q_snow_glacier[:min_len], q_glacier_total[:min_len],
        q_ice_raw[:min_len], q_snow_glacier_raw[:min_len],
        skip_gm=args.skip_gm,
    )

    # 2. Save complete metadata as JSON (all info needed to reproduce)
    metadata = {
        'run_id': RUN_ID,
        'run_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'elapsed_seconds': round(elapsed_total, 1),
        'elapsed_hours': round(elapsed_total / 3600, 2),
        'evaluations': int(nfev),
        
        # Data sources
        'data_sources': {
            'prec_source': args.prec_source,
            'prec_dir': PREC_DIR,
            'temp_dir': TEMP_DIR,
            'evap_dir': EVAP_DIR,
            'glacier_mode': args.glacier_mode,
            'glacier_mask': GLACIER_MASK_PATH,
            'obs_file': OBS_FILE,
        },
        
        # Time configuration
        'time_config': {
            'warmup_start': WARMUP_START,
            'warmup_end': WARMUP_END,
            'calib_start': CALIB_START,
            'calib_end': CALIB_END,
            'valid_start': VALID_START,
            'valid_end': VALID_END,
            'warmup_days': WARMUP_DAYS,
        },
        
        # Basin info
        'basin_info': {
            'catchment_area_km2': CATCHMENT_AREA,
            'valid_cells': len(VALID_CELLS),
            'glacier_cells': int(GLACIER_MASK.sum()),
            'zone_mid_cells': int(ZONE_MID.sum()),
            'zone_high_cells': int(ZONE_HIGH.sum()),
            'area_coef': round(AREA_COEF, 6),
        },
        
        # Optimization settings
        'optimization': {
            'maxiter': args.maxiter,
            'popsize': args.popsize,
            'seed': args.seed,
            'tol': args.tol,
            'workers': args.workers,
            'polish': args.polish,
            'selected_stage': selected_stage,
            'nfev_global': int(nfev_global),
            'nfev_refine': int(nfev_refine),
            'nfev_total': int(nfev),
            'refine': {
                'enabled': bool(do_refine),
                'maxiter': int(REFINE_MAXITER),
                'popsize': int(REFINE_POPSIZE),
                'tol': float(REFINE_TOL),
                'shrink': float(REFINE_SHRINK),
                'min_width': float(REFINE_MIN_WIDTH),
                'seed': int(args.seed + 1),
            },
        },

        # Observation processing
        'obs_processing': {
            'obs_mode_override': OBS_MODE_OVERRIDE,
            'obs_eval_desc': OBS_EVAL_DESC,
            'rule': (
                '若识别到年鉴径流序列存在缺测/推算期（通常 11-4 月为缺测，日值在同一年内呈常数，且 5-10 月占比接近固定值），'
                '则在 NSE/logNSE/PBIAS 等指标计算与目标函数中仅使用 5-10 月（实测期）日值；'
                '其余月份不参与 NSE/logNSE/PBIAS 拟合，仅用于输出与物理诊断；'
                '同时启用弱季节性约束（5-10 月占比下限），避免参数将水量不合理地转移到冬季；'
                '可选启用推算期月体积约束（11-4 月）与冬季冰源约束（12-2 月），用于进一步抑制冬季非物理高流量/融冰。'
            ),
            'obs_used_months': [5, 6, 7, 8, 9, 10] if (OBS_EVAL_MASK is not None and np.any(OBS_OBJ_MASK)) else None,
            'n_obs_used_calib': int(np.sum(np.isfinite(Q_OBS_CALIB))) if Q_OBS_CALIB is not None else None,
            'n_obs_used_valid': int(np.sum(np.isfinite(Q_OBS_VALID))) if Q_OBS_VALID is not None else None,
            'comparison_metrics': 'full_period',
        },

        # Objective configuration (interpretable multi-metric)
        'objective_config': {
            'nse': 'NSE(Q)',
            'log_nse': 'NSE(log1p(Q))',
            'w_val': OBJ_W_VAL,
            'w_log': OBJ_W_LOG,
            'w_pbias_cal': OBJ_W_PBIAS_CAL,
            'w_pbias_val': OBJ_W_PBIAS_VAL,
            'swe_ice_threshold_mm': SWE_ICE_THRESHOLD_MM,
            'muskingum_safe_domain': {
                'dt_day': float(MUSK_DT),
                'x_fixed': float(X_MUSK_FIXED),
                'k_safe_min': float(K_MUSK_MIN_SAFE),
                'k_safe_max': float(K_MUSK_MAX_SAFE),
                'desc': '限制 K_MUSK 使 Muskingum 系数非负，避免出现负系数导致分量路由非线性/不闭合。'
            },
            'seasonality_prior': {
                'enabled': bool(YEARBOOK_MODE and (OBJ_W_SEASON > 0)),
                'mayoct_frac_target': float(YEARBOOK_MAYOCT_FRAC) if (YEARBOOK_MAYOCT_FRAC is not None and np.isfinite(YEARBOOK_MAYOCT_FRAC)) else YEARBOOK_MAYOCT_FRAC_DEFAULT,
                'mayoct_frac_slack': float(YEARBOOK_MAYOCT_FRAC_SLACK),
                'w_season': float(OBJ_W_SEASON),
                'desc': '仅在识别到年鉴缺测/推算期时启用：对模拟 5-10 月占比设置下限(soft penalty)，抑制冬季高流量。'
            },
            'offseason_monthly_prior': {
                'enabled': bool(YEARBOOK_MODE and (OBJ_W_OFFSEASON_MONTHLY > 0)),
                'months': list(OFFSEASON_MONTHS),
                'zero_frac_slack_each_month': float(OFFSEASON_ZERO_FRAC_SLACK),
                'w_offseason_monthly': float(OBJ_W_OFFSEASON_MONTHLY),
                'desc': '仅在年鉴模式启用：利用推算期提供的月尺度体积信息，抑制冬季高流量并约束 4 月回升。'
            },
            'winter_ice_prior': {
                'enabled': bool(YEARBOOK_MODE and (OBJ_W_WINTER_ICE > 0)),
                'winter_months': [12, 1, 2],
                'max_mean_m3s': float(WINTER_ICE_MAX_MEAN_M3S),
                'w_winter_ice': float(OBJ_W_WINTER_ICE),
                'desc': '仅在年鉴模式启用：约束冬季冰源出流接近 0，避免非物理冬季融冰。'
            },
            'gm_volume_prior': {
                'enabled': bool((OBJ_W_GM_VOL > 0) and (GM_SERIES_SIM is not None)),
                'melt_months': [4, 5, 6, 7, 8, 9],
                'ratio_min': float(GM_VOL_RATIO_MIN),
                'ratio_max': float(GM_VOL_RATIO_MAX),
                'w_gm_vol': float(OBJ_W_GM_VOL),
                'desc': '利用 GM_*.tif 的独立冰川融水序列，对融化季(4-9月)模拟冰源体积量级做软约束，防止冰川贡献失真。'
            },
        },
        
        # Optimized parameters
        'optimized_params': {name: round(float(val), 6) for name, val in zip(param_names, result.x)},
        
        # Fixed parameters
        'fixed_params': {k: float(v) for k, v in FIXED.items()},
        
        # Parameter bounds
        'param_bounds': {name: list(bounds) for name, bounds in zip(param_names, PARAM_BOUNDS)},
        
        # Initial states
        'initial_states': {
            'sp': INIT_ST[0], 'sm': INIT_ST[1], 'uz': INIT_ST[2], 
            'lz': INIT_ST[3], 'wc': INIT_ST[4]
        },
        
        # Metrics
        'metrics': {
            'calibration': {
                'nse': round(nse_calib, 4),
                'nse_raw': round(nse_calib_raw, 4),
                'nse_full_period': round(nse_calib_raw, 4),
                'pbias_pct': round(pbias_calib, 2),
                'pbias_raw_pct': round(pbias_calib_raw, 2),
                'pbias_full_period_pct': round(pbias_calib_raw, 2),
                'rmse_m3s': round(rmse_calib, 2),
                'rmse_raw_m3s': round(rmse_calib_raw, 2),
                'rmse_full_period_m3s': round(rmse_calib_raw, 2),
                'rain_fraction': round(rain_frac_calib, 4),
                'snow_noglacier_fraction': round(snow_nogl_frac_calib, 4),
                'glacier_fraction': round(glacier_frac_calib, 4),
                'ice_fraction': round(ice_frac_calib, 4),
                'snow_glacier_fraction': round(snow_glacier_frac_calib, 4),
            },
            'validation': {
                'nse': round(nse_valid, 4),
                'nse_raw': round(nse_valid_raw, 4),
                'nse_full_period': round(nse_valid_raw, 4),
                'pbias_pct': round(pbias_valid, 2),
                'pbias_raw_pct': round(pbias_valid_raw, 2),
                'pbias_full_period_pct': round(pbias_valid_raw, 2),
                'rmse_m3s': round(rmse_valid, 2),
                'rmse_raw_m3s': round(rmse_valid_raw, 2),
                'rmse_full_period_m3s': round(rmse_valid_raw, 2),
                'rain_fraction': round(rain_frac_valid, 4),
                'snow_noglacier_fraction': round(snow_nogl_frac_valid, 4),
                'glacier_fraction': round(glacier_frac_valid, 4),
                'ice_fraction': round(ice_frac_valid, 4),
                'snow_glacier_fraction': round(snow_glacier_frac_valid, 4),
            },
        },
        'metrics_objective_window': {
            'calibration': {
                'nse': round(nse_calib, 4),
                'pbias_pct': round(pbias_calib, 2),
                'rmse_m3s': round(rmse_calib, 2),
            },
            'validation': {
                'nse': round(nse_valid, 4),
                'pbias_pct': round(pbias_valid, 2),
                'rmse_m3s': round(rmse_valid, 2),
            },
        },
        'metrics_full_period': {
            'calibration': {
                'nse': round(nse_calib_raw, 4),
                'pbias_pct': round(pbias_calib_raw, 2),
                'rmse_m3s': round(rmse_calib_raw, 2),
            },
            'validation': {
                'nse': round(nse_valid_raw, 4),
                'pbias_pct': round(pbias_valid_raw, 2),
                'rmse_m3s': round(rmse_valid_raw, 2),
            },
        },

        # Diagnostics (hard-checks + GM comparison; see also diagnostics.json)
        'diagnostics': diag,
    }
    
    metadata_file = os.path.join(run_dir, 'metadata.json')
    with open(metadata_file, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    # 3. Save parameters as simple text file (easy to read)
    param_file = os.path.join(run_dir, 'parameters.txt')
    with open(param_file, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write(f"HBV V6 Calibration Results ({len(param_names)} parameters)\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"Run ID: {RUN_ID}\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Elapsed: {elapsed_total/3600:.2f} hours ({elapsed_total:.0f} seconds)\n")
        f.write(f"Evaluations: {nfev}\n\n")
        
        f.write("-" * 60 + "\n")
        f.write("DATA SOURCES\n")
        f.write("-" * 60 + "\n")
        f.write(f"Precipitation: {args.prec_source} ({PREC_DIR})\n")
        f.write(f"Glacier mode: {args.glacier_mode}\n\n")
        
        f.write("-" * 60 + "\n")
        f.write(f"OPTIMIZED PARAMETERS ({len(param_names)})\n")
        f.write("-" * 60 + "\n")
        for name, val in zip(param_names, result.x):
            f.write(f"{name:15s} = {val:.6f}\n")
        
        f.write("\n" + "-" * 60 + "\n")
        f.write("FIXED PARAMETERS\n")
        f.write("-" * 60 + "\n")
        for name, val in FIXED.items():
            f.write(f"{name:15s} = {val:.6f}\n")
        
        f.write("\n" + "-" * 60 + "\n")
        f.write("PERFORMANCE METRICS\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'Period':<12} {'NSE':>8} {'PBIAS%':>10} {'RMSE':>10} {'Rain%':>8} {'SnowNG%':>9} {'Glac%':>8}\n")
        f.write(f"{'Calibration':<12} {nse_calib:>8.4f} {pbias_calib:>10.2f} {rmse_calib:>10.2f} {rain_frac_calib*100:>8.1f} {snow_nogl_frac_calib*100:>9.1f} {glacier_frac_calib*100:>8.1f}\n")
        f.write(f"{'Validation':<12} {nse_valid:>8.4f} {pbias_valid:>10.2f} {rmse_valid:>10.2f} {rain_frac_valid*100:>8.1f} {snow_nogl_frac_valid*100:>9.1f} {glacier_frac_valid*100:>8.1f}\n")
        f.write("\nNOTE\n")
        f.write("-" * 60 + "\n")
        f.write("SnowNG% = 非冰川区融雪产流占比；Glac% = 冰川区总产流占比(裸冰融 + 冰川区融雪)。\n")
        f.write(f"Calibration: Ice%={ice_frac_calib*100:.1f}%, Glacier-snow%={snow_glacier_frac_calib*100:.1f}%\n")
        f.write(f"Validation : Ice%={ice_frac_valid*100:.1f}%, Glacier-snow%={snow_glacier_frac_valid*100:.1f}%\n")
        if OBS_EVAL_DESC and OBS_EVAL_DESC.startswith("mayoct_only_yearbook"):
            f.write("obs处理: 识别到年鉴缺测/推算期（通常 11-4 月日值在同一年内为常数），\n")
            f.write("         NSE/logNSE/PBIAS 等指标与率定目标仅使用 5-10 月日值（q_obs_obj 其余月份置 NaN）。\n")
            target_frac = YEARBOOK_MAYOCT_FRAC_DEFAULT
            if YEARBOOK_MAYOCT_FRAC is not None and np.isfinite(YEARBOOK_MAYOCT_FRAC):
                target_frac = float(YEARBOOK_MAYOCT_FRAC)
            min_frac = max(0.0, target_frac - YEARBOOK_MAYOCT_FRAC_SLACK)
            f.write(
                f"         弱季节性约束: 要求模拟 5-10 月占比不低于 {min_frac:.3f} "
                f"(target={target_frac:.3f}, slack={YEARBOOK_MAYOCT_FRAC_SLACK:.3f})。\n"
            )
            if OBJ_W_OFFSEASON_MONTHLY > 0:
                f.write(
                    f"         推算期月体积约束: 11-4 月按观测月体积/占比做软惩罚 "
                    f"(w={OBJ_W_OFFSEASON_MONTHLY:.1f}, zero_slack={OFFSEASON_ZERO_FRAC_SLACK:.3f})。\n"
                )
            if OBJ_W_WINTER_ICE > 0:
                f.write(
                    f"         冬季冰源约束: 12-2 月冰源出流均值不高于 {WINTER_ICE_MAX_MEAN_M3S:.2f} m3/s "
                    f"(w={OBJ_W_WINTER_ICE:.1f})。\n"
                )
        else:
            f.write("obs处理: 默认使用全年的日流量参与 NSE/logNSE/PBIAS 指标与率定目标。\n")

        if (OBJ_W_GM_VOL > 0) and (GM_SERIES_SIM is not None):
            f.write(
                f"GM约束: 融化季(4-9月)模拟冰源体积/GM体积应位于 "
                f"[{GM_VOL_RATIO_MIN:.2f}, {GM_VOL_RATIO_MAX:.2f}] (w={OBJ_W_GM_VOL:.1f})。\n"
            )
        f.write(f"Full-period reference: 率定期 NSE={nse_calib_raw:.4f}, PBIAS={pbias_calib_raw:+.2f}%, RMSE={rmse_calib_raw:.2f}\n")
        f.write(f"Full-period reference: 验证期 NSE={nse_valid_raw:.4f}, PBIAS={pbias_valid_raw:+.2f}%, RMSE={rmse_valid_raw:.2f}\n")

    # 4. Copy log file to run directory
    if os.path.exists(LOG_FILE):
        import shutil
        shutil.copy(LOG_FILE, os.path.join(run_dir, 'optimization.log'))
    if PROGRESS_FILE_MAIN and os.path.exists(PROGRESS_FILE_MAIN):
        import shutil
        shutil.copy(PROGRESS_FILE_MAIN, os.path.join(run_dir, 'progress.csv'))
    if PROGRESS_FILE_REFINE and os.path.exists(PROGRESS_FILE_REFINE):
        import shutil
        shutil.copy(PROGRESS_FILE_REFINE, os.path.join(run_dir, 'progress_refine.csv'))

    print(f"\n" + "=" * 60)
    print(f"RESULTS SAVED TO: {run_dir}")
    print("=" * 60)
    print(f"  - simulation.csv    : 模拟结果时间序列")
    print(f"  - simulation_detailed.csv : 分量细分(含冰川雪融/裸冰融)与原始(未汇流)序列")
    print(f"  - metadata.json     : 完整元数据 (可复现)")
    print(f"  - diagnostics.json  : 物理自检 + GM 对照 (硬约束检查)")
    if os.path.exists(os.path.join(run_dir, 'cryo_diagnostics.csv')):
        print(f"  - cryo_diagnostics.csv : GM 与冰川分量对照时序")
    print(f"  - parameters.txt    : 参数汇总 (易读)")
    print(f"  - optimization.log  : 优化日志")
    print(f"  - progress.csv      : 优化进度")
    if os.path.exists(os.path.join(run_dir, 'progress_refine.csv')):
        print(f"  - progress_refine.csv : 二次精修优化进度")


if __name__ == "__main__":
    main()

