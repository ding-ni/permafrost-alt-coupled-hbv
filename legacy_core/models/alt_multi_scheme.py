# -*- coding: utf-8 -*-
"""
实验2多方案对比：多年冻土-HBV耦合模型

基于物理机制设计多种ALT参数修正方案，一次运行对比所有方案效果。

物理机制分析：
1. ALT (活动层厚度) 决定了季节性融化土层的深度
2. ALT越厚 → 土壤蓄水空间越大 → FC应该增大（正向关系）
3. ALT越厚 → 深层排水通道越长 → K2应该增大（正向关系）
4. ALT越厚 → 土壤不易饱和 → K0（快速径流）应该减小（反向关系）

修正方案设计：
- 方案0: V6基线（无ALT修正）
- 方案1: 原始方案 (s=ALT_ref/ALT, FC/PERC×s)
- 方案2: 原始+放大 (s=ALT_ref/ALT, α=2)
- 方案3: 正向修正 (s=ALT/ALT_ref, FC/PERC×s) ← 物理上更合理
- 方案4: 正向+放大 (s=ALT/ALT_ref, α=2)
- 方案5: 混合方案 (FC/K2正向, K0反向)
- 方案6: 全参数修正 (FC/PERC/K2/K0/K)
"""
import sys
import os
import time
import numpy as np
import pandas as pd
import json
from datetime import datetime
from glob import glob
from math import sin, radians
import rasterio
from rasterio.warp import reproject, Resampling
from numba import njit
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 配置
# =============================================================================
PROJECT_ROOT = r'.'
PREC_DIR = os.path.join(PROJECT_ROOT, 'data', 'aligned_masked', 'prec')
TEMP_DIR = os.path.join(PROJECT_ROOT, 'data', 'aligned_masked', 'temp')
EVAP_DIR = os.path.join(PROJECT_ROOT, 'data', 'aligned_masked', 'evap')
OBS_FILE = os.path.join(PROJECT_ROOT, 'data', 'observed', 'discharge_example.csv')
FLOW_ACC_PATH = os.path.join(PROJECT_ROOT, 'data', 'gis', 'flow_accumulation_masked.tif')
GLACIER_MASK_PATH = os.path.join(PROJECT_ROOT, 'data', 'gis', 'glacier_mask.tif')
ALT_DIR = os.path.join(PROJECT_ROOT, 'data', 'alt', 'basin')
SF_DIR = os.path.join(PROJECT_ROOT, 'data', 'alt', 'scale_factors')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results', 'runs')

# ALT_ref 默认文件名（可在 CLI 通过 --alt-ref 覆盖）
# - 为避免验证期信息泄露，若存在 ALT_ref_2006_2017.tif 则优先使用
ALT_REF_FILE = "ALT_ref.tif"

CATCHMENT_AREA = 15924
SWE_ICE_THRESHOLD_MM = 10.0
EPS = 1e-12

WARMUP_START = '2006-01-01'
WARMUP_END = '2008-12-31'
CALIB_START = '2009-01-01'
CALIB_END = '2017-12-31'
VALID_START = '2018-01-01'
VALID_END = '2020-12-31'
SIM_START = CALIB_START
SIM_END = VALID_END
OBS_MODE_OVERRIDE = 'auto'  # auto | full_year | may_oct

MUSK_DT = 1.0
X_MUSK_FIXED = 0.2

# =============================================================================
# V6 率定参数
# =============================================================================
V6_PARAMS = {
    'TT': -1.854508,
    'FC': 1401.862945,
    'BETA': 1.250633,
    'LP': 0.824162,
    'RFCF': 0.856409,
    'SFCF': 0.506925,
    'CFR': 0.066875,
    'CWH': 0.036316,
    'CFMAX_mid': 4.988928,
    'CFMAX_high': 3.875559,
    'K': 0.073742,
    'K0': 0.073742,  # 使用K作为K0的基准（V6中K≈K0）
    'K1': 0.051772,
    'K2': 0.001354,
    'PERC': 0.061616,
    'K_MUSK': 0.984463,
    'ICE_FACTOR': 1.388001,
}

FIXED = {
    'ALPHA': 0.0,
    'UZL': 30.0,
    'X_MUSK': X_MUSK_FIXED,
}

# 参数物理约束
PARAM_BOUNDS = {
    'FC': (200.0, 2500.0),
    'PERC': (0.01, 0.5),
    'K2': (0.0002, 0.02),
    'K0': (0.02, 0.5),
    'K': (0.02, 0.5),
}

INIT_ST = np.array([0.0, 5.0, 0.0, 0.0, 0.0], dtype=np.float64)

# =============================================================================
# 方案定义（基于实验结果优化）
# =============================================================================
# 实验发现：
# 1. inverse公式(s=ALT_ref/ALT)比forward公式效果好
# 2. 放大因子α=2有效
# 3. forward方案(S3-S6)反而不如基线，已删除
#
# 进一步探索：
# - 不同放大因子 α = 1.5, 2.0, 2.5, 3.0
# - 不同参数组合

SCHEMES = {
    'baseline': {
        'description': 'V6基线（无ALT修正）',
        'use_alt': False,
    },
    'S1_fc_perc': {
        'description': 's=ALT_ref/ALT, α=1, FC+PERC',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 1.0,
        'params': ['FC', 'PERC'],
    },
    'S2_fc_perc_k2_a1': {
        'description': 's=ALT_ref/ALT, α=1, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 1.0,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S3_fc_perc_k2_a15': {
        'description': 's=ALT_ref/ALT, α=1.5, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 1.5,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S4_fc_perc_k2_a2': {
        'description': 's=ALT_ref/ALT, α=2, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 2.0,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S5_fc_perc_k2_a25': {
        'description': 's=ALT_ref/ALT, α=2.5, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 2.5,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S6_fc_perc_k2_a3': {
        'description': 's=ALT_ref/ALT, α=3, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 3.0,
        'params': ['FC', 'PERC', 'K2'],
    },
    'fc_only_alpha2': {
        'description': 's=ALT_ref/ALT, α=2, 仅FC',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 2.0,
        'params': ['FC'],
    },
    'k2_only_alpha2': {
        'description': 's=ALT_ref/ALT, α=2, 仅K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 2.0,
        'params': ['K2'],
    },
    'S9_fc_perc_k2_a4': {
        'description': 's=ALT_ref/ALT, α=4, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 4.0,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S10_fc_perc_k2_a5': {
        'description': 's=ALT_ref/ALT, α=5, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 5.0,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S11_fc_perc_k2_a6': {
        'description': 's=ALT_ref/ALT, α=6, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 6.0,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S12_fc_perc_k2_a8': {
        'description': 's=ALT_ref/ALT, α=8, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 8.0,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S13_fc_perc_k2_a10': {
        'description': 's=ALT_ref/ALT, α=10, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 10.0,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S21_fc_perc_k2_a10_clip': {
        'description': 's=ALT_ref/ALT, α=10, FC+PERC+K2, factor_clip[0.2,3.0]',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 10.0,
        'factor_clip': (0.2, 3.0),
        'params': ['FC', 'PERC', 'K2'],
    },
    'S22_fc_perc_k2_a12_clip': {
        'description': 's=ALT_ref/ALT, α=12, FC+PERC+K2, factor_clip[0.2,3.0]',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 12.0,
        'factor_clip': (0.2, 3.0),
        'params': ['FC', 'PERC', 'K2'],
    },
    'S23_fc_perc_k2_a14_clip': {
        'description': 's=ALT_ref/ALT, α=14, FC+PERC+K2, factor_clip[0.2,3.0]',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 14.0,
        'factor_clip': (0.2, 3.0),
        'params': ['FC', 'PERC', 'K2'],
    },
    'S24_fc_perc_k2_exp4': {
        'description': 's=ALT_ref/ALT, exp=4, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'amp_mode': 'exp',
        'alpha': 4.0,
        'params': ['FC', 'PERC', 'K2'],
    },
    'fc_perc_k2_alpha10': {
        'description': 'post-2018 only: s=ALT_ref/ALT, α=10, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 10.0,
        'start_year': 2018,
        'params': ['FC', 'PERC', 'K2'],
    },
    'fc_perc_k2_alpha6': {
        'description': 'post-2018 only: s=ALT_ref/ALT, α=6, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 6.0,
        'start_year': 2018,
        'params': ['FC', 'PERC', 'K2'],
    },
    'fc_perc_k2_alpha8': {
        'description': 'post-2018 only: s=ALT_ref/ALT, α=8, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 8.0,
        'start_year': 2018,
        'params': ['FC', 'PERC', 'K2'],
    },
    'fc_perc_k2_alpha12': {
        'description': 'post-2018 only: s=ALT_ref/ALT, α=12, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 12.0,
        'start_year': 2018,
        'params': ['FC', 'PERC', 'K2'],
    },
    'fc_perc_k2_alpha14': {
        'description': 'post-2018 only: s=ALT_ref/ALT, α=14, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 14.0,
        'start_year': 2018,
        'params': ['FC', 'PERC', 'K2'],
    },
    'fc_perc_k2_alpha10_clipped': {
        'description': 'post-2018 only: s=ALT_ref/ALT, α=10, FC+PERC+K2, factor_clip[0.2,3.0]',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 10.0,
        'start_year': 2018,
        'factor_clip': (0.2, 3.0),
        'params': ['FC', 'PERC', 'K2'],
    },

    # ---------------------------------------------------------------------
    # More defensible (physics-consistent sign) post-2018 families:
    # Aim: ALT deepening -> FC/PERC/K2 increase. Implemented in two ways:
    # 1) inverse scale factor but invert params (divide by factor) to flip the sign.
    # 2) forward scale factor (ALT/ALT_ref) with positive coupling.
    # ---------------------------------------------------------------------
    'physically_inverted_fc_perc_k2_alpha6': {
        'description': 'post-2018 only (phys-sign): s=ALT_ref/ALT, α=6, FC+PERC+K2 (inverse_params)',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 6.0,
        'start_year': 2018,
        'factor_clip': (0.2, 3.0),
        'inverse_params': ['FC', 'PERC', 'K2'],
        'params': ['FC', 'PERC', 'K2'],
    },
    'physically_inverted_fc_perc_k2_alpha10': {
        'description': 'post-2018 only (phys-sign): s=ALT_ref/ALT, α=10, FC+PERC+K2 (inverse_params)',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 10.0,
        'start_year': 2018,
        'factor_clip': (0.2, 3.0),
        'inverse_params': ['FC', 'PERC', 'K2'],
        'params': ['FC', 'PERC', 'K2'],
    },
    'forward_fc_perc_k2_alpha3': {
        'description': 'post-2018 only (forward): s=ALT/ALT_ref, α=3, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'forward',
        'alpha': 3.0,
        'start_year': 2018,
        'factor_clip': (0.2, 3.0),
        'params': ['FC', 'PERC', 'K2'],
    },
    'forward_fc_perc_k2_exponent2': {
        'description': 'post-2018 only (forward-exp): s=ALT/ALT_ref, exp=2, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'forward',
        'amp_mode': 'exp',
        'alpha': 2.0,
        'start_year': 2018,
        'params': ['FC', 'PERC', 'K2'],
    },
    'S28_ramp2015_2020_fc_perc_k2_a10': {
        'description': 'ramp 2015-2020: s=ALT_ref/ALT, α=10, FC+PERC+K2',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 10.0,
        'year_ramp': (2015, 2020),
        'params': ['FC', 'PERC', 'K2'],
    },
    'S25_zoned_high_only_a10': {
        'description': 's=ALT_ref/ALT, α=10, FC+PERC+K2, only high zone',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 10.0,
        'params': ['FC', 'PERC', 'K2'],
        'zone_weights': {'high': 1.0, 'mid': 0.0, 'low': 0.0},
    },
    'S26_zoned_high15_mid10_low05_a10': {
        'description': 's=ALT_ref/ALT, α=10, FC+PERC+K2, zoned weights (1.5/1.0/0.5)',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 10.0,
        'params': ['FC', 'PERC', 'K2'],
        'zone_weights': {'high': 1.5, 'mid': 1.0, 'low': 0.5},
    },

    # ---------------------------------------------------------------------
    # Extended inverse-family (empirically best formula) with quick-response K:
    # K controls fast outflow above UZL (HBV K0-like). Adding K may improve peaks/timing
    # without needing full re-calibration.
    # ---------------------------------------------------------------------
    'S17_inv_fc_perc_k2_k_a5': {
        'description': 's=ALT_ref/ALT, α=5, FC+PERC+K2+K',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 5.0,
        'params': ['FC', 'PERC', 'K2', 'K'],
    },
    'S18_inv_fc_perc_k2_k_inv_a5': {
        'description': 's=ALT_ref/ALT, α=5, FC+PERC+K2, K inverse',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 5.0,
        'params': ['FC', 'PERC', 'K2', 'K'],
        'inverse_params': ['K'],
    },
    'S19_inv_fc_perc_k2_k_a6': {
        'description': 's=ALT_ref/ALT, α=6, FC+PERC+K2+K',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': 6.0,
        'params': ['FC', 'PERC', 'K2', 'K'],
    },
    'S20_inv_weighted_fcdom': {
        'description': 's=ALT_ref/ALT, α={FC:8,PERC:5,K2:2,K:3}, FC-dominant',
        'use_alt': True,
        'scale_formula': 'inverse',
        'alpha': {'FC': 8.0, 'PERC': 5.0, 'K2': 2.0, 'K': 3.0, 'default': 1.0},
        'params': ['FC', 'PERC', 'K2', 'K'],
    },

    # ---------------------------------------------------------------------
    # Physically interpretable family:
    # ALT deepening (thaw) -> higher subsurface storage/percolation (FC/PERC/K2 up)
    # ALT shallowing (stronger frost barrier) -> stronger quick response (K up)
    # Implemented as forward scale factor for FC/PERC/K2 and inverse for K.
    # ---------------------------------------------------------------------
    'S14_phys_forward_fc_perc_k_inv_exp2': {
        'description': 'phys: s=ALT/ALT_ref, exp=2, FC+PERC, K inverse',
        'use_alt': True,
        'scale_formula': 'forward',
        'amp_mode': 'exp',
        'alpha': 2.0,
        'params': ['FC', 'PERC', 'K'],
        'inverse_params': ['K'],
    },
    'S15_phys_forward_fc_perc_k2_k_inv_exp2': {
        'description': 'phys: s=ALT/ALT_ref, exp=2, FC+PERC+K2, K inverse',
        'use_alt': True,
        'scale_formula': 'forward',
        'amp_mode': 'exp',
        'alpha': 2.0,
        'params': ['FC', 'PERC', 'K2', 'K'],
        'inverse_params': ['K'],
    },
    'S16_phys_forward_weighted_exp': {
        'description': 'phys-weighted: s=ALT/ALT_ref, exp={FC:2.5,PERC:2,K2:1.5,K:2.5}, K inverse',
        'use_alt': True,
        'scale_formula': 'forward',
        'amp_mode': 'exp',
        'alpha': {'FC': 2.5, 'PERC': 2.0, 'K2': 1.5, 'K': 2.5, 'default': 1.0},
        'params': ['FC', 'PERC', 'K2', 'K'],
        'inverse_params': ['K'],
    },
}

# =============================================================================
# 全局变量
# =============================================================================
PREC_3D = None
TEMP_3D = None
ET_3D = None
LL_TEMP_3D = None
VALID_CELLS = None
CELL_SCALE = None
GLACIER_MASK = None
Q_OBS_FULL = None
Q_OBS_OBJ_FULL = None
WARMUP_DAYS = None
SIM_DATES = None
ZONE_MID = None
ZONE_HIGH = None
OBS_EVAL_MASK = None
OBS_EVAL_DESC = None

# ALT数据
ALT_REF = None  # 多年平均
ALT_YEARLY = {}  # {year: alt_grid}


# =============================================================================
# Numba函数
# =============================================================================
@njit(cache=True, fastmath=True, nogil=True)
def hbv_cell_full(prec, temp, et, ll_temp, par, init_st, is_glacier, ice_factor):
    """
    单像元HBV模拟（支持所有参数空间变化）

    par数组顺序:
    [0]TT, [1]RFCF, [2]SFCF, [3]CFMAX, [4]CWH, [5]CFR,
    [6]FC, [7]BETA, [8]ALPHA, [9]LP, [10]K, [11]K1, [12]K2, [13]UZL, [14]PERC
    """
    n = len(prec)
    q_uz = np.zeros(n, dtype=np.float32)
    q_lz = np.zeros(n, dtype=np.float32)
    q_rain = np.zeros(n, dtype=np.float32)
    q_snow = np.zeros(n, dtype=np.float32)
    q_ice = np.zeros(n, dtype=np.float32)
    final_st = np.zeros(5, dtype=np.float64)

    tt, rfcf, sfcf, cfmax, cwh, cfr = par[0], par[1], par[2], par[3], par[4], par[5]
    fc, beta, e_corr, lp = par[6], par[7], par[8], par[9]
    k, k1, k2, uzl, perc = par[10], par[11], par[12], par[13], par[14]

    fc = max(fc, 10.0)
    beta = max(min(beta, 10.0), 0.1)
    lp = max(min(lp, 0.99), 0.01)

    sp, sm, uz, lz, wc = init_st[0], init_st[1], init_st[2], init_st[3], init_st[4]
    # 关键修复：将传递过来的存储量赋给rain组分，保证水量守恒
    uz_r, uz_s, uz_i = uz, 0.0, 0.0
    lz_r, lz_s, lz_i = lz, 0.0, 0.0

    for i in range(n):
        p = prec[i]
        t = temp[i]
        e = et[i]
        tm = ll_temp[i]

        if np.isnan(p) or np.isnan(t) or np.isnan(e):
            continue

        p = max(p, 0.0)
        e = max(e, 0.0)

        if t <= tt:
            rf, sf = 0.0, p * sfcf
        else:
            rf, sf = p * rfcf, 0.0

        melt = 0.0
        ice_melt = 0.0
        if t > tt:
            ddt = t - tt
            avail_snow = sp + sf
            melt_pot_snow = cfmax * ddt
            melt = min(melt_pot_snow, avail_snow)
            sp = max(avail_snow - melt, 0.0)

            if is_glacier and (sp <= SWE_ICE_THRESHOLD_MM) and (melt_pot_snow > 0.0):
                f_snow = melt / (melt_pot_snow + EPS)
                f_snow = max(0.0, min(1.0, f_snow))
                melt_pot_ice = (cfmax * ice_factor) * ddt
                ice_melt = max(melt_pot_ice * (1.0 - f_snow), 0.0)

            wc_int = wc + melt + rf + ice_melt
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

        input_total = rf + melt + ice_melt
        den_in = input_total + EPS
        frac_r = rf / den_in
        frac_s = melt / den_in
        frac_i = ice_melt / den_in

        sm_ratio = min(max(sm / fc, 0.0), 1.0)
        r = (sm_ratio ** beta) * inf
        r_r, r_s, r_i = r * frac_r, r * frac_s, r * frac_i

        ep_adj = max((1.0 + (t - tm) * e_corr) * e, 0.0)
        lp_fc = lp * fc
        ea = min(ep_adj, (sm / lp_fc) * ep_adj) if lp_fc > 0.001 else ep_adj
        ea = min(ea, sm)

        uz_r += r_r
        uz_s += r_s
        uz_i += r_i
        uz_int = uz_r + uz_s + uz_i
        sm = max(min(sm + inf - r - ea, fc), 0.0)

        perc_actual = min(perc, uz_int)
        den_uz = uz_int + EPS
        frac_ur, frac_us, frac_ui = uz_r / den_uz, uz_s / den_uz, uz_i / den_uz
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
            q0, q1 = uz_int2 * 0.67, uz_int2 * 0.33

        den_uz2 = uz_int2 + EPS
        frac_ur2, frac_us2, frac_ui2 = uz_r / den_uz2, uz_s / den_uz2, uz_i / den_uz2

        q0_r, q0_s, q0_i = q0 * frac_ur2, q0 * frac_us2, q0 * frac_ui2
        q1_r, q1_s, q1_i = q1 * frac_ur2, q1 * frac_us2, q1 * frac_ui2

        uz_r = max(uz_r - (q0_r + q1_r), 0.0)
        uz_s = max(uz_s - (q0_s + q1_s), 0.0)
        uz_i = max(uz_i - (q0_i + q1_i), 0.0)

        lz_int = lz_r + lz_s + lz_i
        q2 = k2 * lz_int
        if q2 > lz_int:
            q2 = lz_int

        den_lz = lz_int + EPS
        frac_lr, frac_ls, frac_li = lz_r / den_lz, lz_s / den_lz, lz_i / den_lz

        q2_r, q2_s, q2_i = q2 * frac_lr, q2 * frac_ls, q2 * frac_li

        lz_r = max(lz_r - q2_r, 0.0)
        lz_s = max(lz_s - q2_s, 0.0)
        lz_i = max(lz_i - q2_i, 0.0)

        # Clip storages
        uz_sum = uz_r + uz_s + uz_i
        if uz_sum > 10000.0:
            scale = 10000.0 / uz_sum
            uz_r, uz_s, uz_i = uz_r * scale, uz_s * scale, uz_i * scale

        lz_sum = lz_r + lz_s + lz_i
        if lz_sum > 10000.0:
            scale = 10000.0 / lz_sum
            lz_r, lz_s, lz_i = lz_r * scale, lz_s * scale, lz_i * scale

        q_uz[i] = q0 + q1
        q_lz[i] = q2
        q_rain[i] = q0_r + q1_r + q2_r
        q_snow[i] = q0_s + q1_s + q2_s
        q_ice[i] = q0_i + q1_i + q2_i

    final_st[0] = sp
    final_st[1] = sm
    final_st[2] = uz_r + uz_s + uz_i
    final_st[3] = lz_r + lz_s + lz_i
    final_st[4] = wc

    return q_uz, q_lz, q_rain, q_snow, q_ice, final_st


@njit(cache=True, fastmath=True, nogil=True)
def run_all_cells_scheme(prec_3d, temp_3d, et_3d, ll_temp_3d, par_base,
                         cfmax_grid, fc_grid, perc_grid, k2_grid, k_grid,
                         init_st_grid, valid_cells, cell_scale, glacier_mask, ice_factor):
    """运行所有像元（支持FC/PERC/K2/K空间变化）"""
    n_cells = valid_cells.shape[0]
    ts = prec_3d.shape[2]

    q_total = np.zeros(ts, dtype=np.float64)
    q_rain = np.zeros(ts, dtype=np.float64)
    q_snow = np.zeros(ts, dtype=np.float64)
    q_ice = np.zeros(ts, dtype=np.float64)
    end_st_grid = np.zeros((n_cells, 5), dtype=np.float64)

    for idx in range(n_cells):
        x, y = valid_cells[idx, 0], valid_cells[idx, 1]
        scale = cell_scale[idx]

        prec = prec_3d[x, y, :]
        temp = temp_3d[x, y, :]
        et = et_3d[x, y, :]
        ll_temp = ll_temp_3d[x, y, :]

        par = par_base.copy()
        par[3] = cfmax_grid[x, y]    # CFMAX
        par[6] = fc_grid[x, y]       # FC
        par[10] = k_grid[x, y]       # K (也用作K0)
        par[12] = k2_grid[x, y]      # K2
        par[14] = perc_grid[x, y]    # PERC

        is_glacier = glacier_mask[x, y]
        init_st = init_st_grid[idx, :]

        q_uz, q_lz, q_r, q_s, q_i, final_st = hbv_cell_full(
            prec, temp, et, ll_temp, par, init_st, is_glacier, ice_factor
        )

        end_st_grid[idx, :] = final_st

        for t in range(ts):
            q_total[t] += (q_uz[t] + q_lz[t]) * scale
            q_rain[t] += q_r[t] * scale
            q_snow[t] += q_s[t] * scale
            q_ice[t] += q_i[t] * scale

    return q_total, q_rain, q_snow, q_ice, end_st_grid


@njit(cache=True, fastmath=True, nogil=True)
def muskingum_route(q_in, k, x):
    n = len(q_in)
    q_out = np.zeros(n, dtype=np.float64)
    q_out[0] = q_in[0]
    dt = 1.0
    denom = 2*k*(1-x) + dt
    c0, c1, c2 = (dt - 2*k*x) / denom, (dt + 2*k*x) / denom, (2*k*(1-x) - dt) / denom
    for i in range(1, n):
        q_out[i] = max(c0 * q_in[i] + c1 * q_in[i-1] + c2 * q_out[i-1], 0.0)
    return q_out


@njit(cache=True, nogil=True)
def nse_calc(obs, sim):
    n = len(obs)
    sum_obs, count = 0.0, 0
    for i in range(n):
        if not (np.isnan(obs[i]) or np.isnan(sim[i])):
            sum_obs += obs[i]
            count += 1
    if count == 0:
        return -999.0
    mean_obs = sum_obs / count
    ss_err, ss_tot = 0.0, 0.0
    for i in range(n):
        if not (np.isnan(obs[i]) or np.isnan(sim[i])):
            ss_err += (obs[i] - sim[i]) ** 2
            ss_tot += (obs[i] - mean_obs) ** 2
    return 1.0 - ss_err / ss_tot if ss_tot > 0 else -999.0


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
    all_files = glob(os.path.join(directory, '*.tif'))
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


def resample_alt_to_ref(alt_file, ref_profile, ref_shape):
    """将ALT文件重采样到参考网格（与ALT_ref一致）"""
    with rasterio.open(alt_file) as src:
        alt_data = src.read(1)
        src_nodata = src.nodata if src.nodata is not None else -9999

        # 将nodata转为特殊值用于重采样
        alt_data = np.where(alt_data == src_nodata, -9999, alt_data)

        # 如果尺寸不匹配，需要重采样
        if alt_data.shape != ref_shape:
            dst_data = np.zeros(ref_shape, dtype=np.float32)

            reproject(
                source=alt_data,
                destination=dst_data,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=ref_profile['transform'],
                dst_crs=ref_profile['crs'],
                resampling=Resampling.bilinear,
                src_nodata=-9999,
                dst_nodata=-9999
            )

            # 转回nan
            alt_data = np.where(dst_data == -9999, np.nan, dst_data)
        else:
            alt_data = np.where(alt_data == -9999, np.nan, alt_data)

    return alt_data.astype(np.float64)


def load_alt_data():
    """加载ALT数据（重采样到统一网格）"""
    global ALT_REF, ALT_YEARLY

    print("\n加载ALT数据...")

    # 加载ALT_ref并获取参考网格信息
    alt_ref_file = ALT_REF_FILE
    if not os.path.isabs(alt_ref_file):
        alt_ref_file = os.path.join(SF_DIR, alt_ref_file)
    ref_profile = None
    ref_shape = None

    if os.path.exists(alt_ref_file):
        with rasterio.open(alt_ref_file) as src:
            ALT_REF = src.read(1).astype(np.float64)
            nodata = src.nodata
            if nodata is not None:
                ALT_REF = np.where(ALT_REF == nodata, np.nan, ALT_REF)
            ref_profile = src.profile.copy()
            ref_shape = (src.height, src.width)
        print(f"   ALT_ref({os.path.basename(alt_ref_file)}): shape={ref_shape}, "
              f"mean={np.nanmean(ALT_REF):.2f}m, range=[{np.nanmin(ALT_REF):.2f}, {np.nanmax(ALT_REF):.2f}]")
    else:
        print(f"   [WARN] ALT_ref 文件不存在: {alt_ref_file}")
        return

    # 加载逐年ALT（重采样到参考网格）
    for year in range(2006, 2021):
        alt_file = os.path.join(ALT_DIR, f"ALT_{year}.tif")
        if os.path.exists(alt_file):
            alt_data = resample_alt_to_ref(alt_file, ref_profile, ref_shape)
            ALT_YEARLY[year] = alt_data

    print(f"   加载了 {len(ALT_YEARLY)} 年ALT数据（已重采样到 {ref_shape}）")


def detect_obs_eval_mask(q_obs_full):
    full_mask = np.ones(len(SIM_DATES), dtype=bool)
    may_oct_mask = (SIM_DATES.month >= 5) & (SIM_DATES.month <= 10)
    mode = str(OBS_MODE_OVERRIDE or 'auto').strip().lower()

    if mode == 'full_year':
        return full_mask, 'full_year_forced'
    if mode == 'may_oct':
        return may_oct_mask, 'mayoct_forced'
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
            return may_oct_mask, f"mayoct_only_yearbook(frac={frac_mean:.3f}±{frac_std:.3g})"
    except Exception as exc:
        print(f"      [WARN] Obs mode detection failed, fallback to full-year eval: {exc}")

    return full_mask, 'full_year'


def load_all_data():
    global PREC_3D, TEMP_3D, ET_3D, LL_TEMP_3D, VALID_CELLS, CELL_SCALE
    global GLACIER_MASK, Q_OBS_FULL, Q_OBS_OBJ_FULL, WARMUP_DAYS, SIM_DATES, ZONE_MID, ZONE_HIGH
    global OBS_EVAL_MASK, OBS_EVAL_DESC

    print("="*70)
    print("加载数据...")
    print("="*70)

    print("[1/6] 加载降水...")
    PREC_3D, transform = load_raster_stack(PREC_DIR, WARMUP_START, SIM_END)
    rows, cols, ts = PREC_3D.shape
    print(f"      Shape: {rows} x {cols} x {ts}")

    print("[2/6] 加载温度...")
    TEMP_3D, _ = load_raster_stack(TEMP_DIR, WARMUP_START, SIM_END)

    print("[3/6] 加载蒸散发...")
    ET_3D, _ = load_raster_stack(EVAP_DIR, WARMUP_START, SIM_END)

    LL_TEMP_3D = np.nanmean(TEMP_3D, axis=2, keepdims=True)
    LL_TEMP_3D = np.broadcast_to(LL_TEMP_3D, TEMP_3D.shape).astype(np.float32)

    PREC_3D = np.nan_to_num(PREC_3D, nan=0.0)
    ET_3D = np.nan_to_num(ET_3D, nan=0.0)
    TEMP_3D = np.where(np.isnan(TEMP_3D), LL_TEMP_3D, TEMP_3D)

    print("[4/6] 加载流量累积...")
    with rasterio.open(FLOW_ACC_PATH) as src:
        flow_acc = src.read(1).astype(np.float32)
        nodata = src.nodata
        if nodata is not None:
            flow_acc[flow_acc == nodata] = np.nan

    valid_mask = ~np.isnan(flow_acc)
    VALID_CELLS = np.argwhere(valid_mask).astype(np.int64)
    print(f"      有效像元: {len(VALID_CELLS)}")

    print("[4.5/6] 加载高程分区...")
    zone_mid_file = os.path.join(PROJECT_ROOT, 'data/gis/elevation_zone_mid.tif')
    zone_high_file = os.path.join(PROJECT_ROOT, 'data/gis/elevation_zone_high.tif')
    with rasterio.open(zone_mid_file) as src:
        ZONE_MID = src.read(1).astype(bool)
    with rasterio.open(zone_high_file) as src:
        ZONE_HIGH = src.read(1).astype(bool)

    # 像元面积
    px_area = np.zeros((rows, cols), dtype=np.float64)
    row_areas = compute_row_areas_km2(transform, rows)
    for i in range(rows):
        px_area[i, :] = row_areas[i]

    px_tot_area = np.sum(px_area[valid_mask])
    area_coef = CATCHMENT_AREA / px_tot_area
    CELL_SCALE = (px_area[VALID_CELLS[:,0], VALID_CELLS[:,1]] * area_coef / 86.4).astype(np.float64)

    print("[5/6] 加载实测流量...")
    SIM_DATES = pd.date_range(start=SIM_START, end=SIM_END, freq='D')
    if os.path.exists(OBS_FILE):
        df_obs = pd.read_csv(OBS_FILE)
        df_obs['date'] = pd.to_datetime(df_obs['date'])
        df_obs.set_index('date', inplace=True)
        Q_OBS_FULL = df_obs.reindex(SIM_DATES)['discharge (m3/s)'].values.astype(np.float64)
        print(f"      实测均值: {np.nanmean(Q_OBS_FULL):.2f} m3/s")
    else:
        Q_OBS_FULL = np.zeros(len(SIM_DATES))

    OBS_EVAL_MASK, OBS_EVAL_DESC = detect_obs_eval_mask(Q_OBS_FULL)
    Q_OBS_OBJ_FULL = Q_OBS_FULL.copy()
    Q_OBS_OBJ_FULL[~OBS_EVAL_MASK] = np.nan
    print(f"      观测口径: {OBS_EVAL_DESC}")
    if OBS_EVAL_DESC.startswith('mayoct'):
        print("      目标窗口 = 5-10月；同时输出全时段参考指标。")
    else:
        print("      目标窗口 = 全年。")

    print("[6/6] 加载冰川掩膜...")
    if os.path.exists(GLACIER_MASK_PATH):
        with rasterio.open(GLACIER_MASK_PATH) as src:
            glacier_raw = src.read(1)
            nodata = src.nodata
            if nodata is not None:
                glacier_raw = np.where(glacier_raw == nodata, 0, glacier_raw)
        GLACIER_MASK = glacier_raw.astype(bool)
        print(f"      冰川像元: {GLACIER_MASK.sum()}")
    else:
        GLACIER_MASK = np.zeros((rows, cols), dtype=bool)

    WARMUP_DAYS = (pd.to_datetime(CALIB_START) - pd.to_datetime(WARMUP_START)).days

    # 加载ALT数据
    load_alt_data()

    print("="*70)


# =============================================================================
# 参数修正函数
# =============================================================================
def compute_scale_factor(year, scheme):
    """计算指定年份的尺度因子"""
    if not scheme.get('use_alt', False):
        return None, None, None

    if year not in ALT_YEARLY or ALT_REF is None:
        return None, None, None

    alt_y = ALT_YEARLY[year]
    formula = scheme.get('scale_formula', 'forward')
    clip_low, clip_high = scheme.get('clip', (0.5, 2.0))

    with np.errstate(divide='ignore', invalid='ignore'):
        valid_mask = (~np.isnan(ALT_REF)) & (~np.isnan(alt_y)) & (ALT_REF > 0.01) & (alt_y > 0.01)
        if formula == 'inverse':
            # 原始公式: s = ALT_ref / ALT_y
            s_raw = np.where(valid_mask, ALT_REF / alt_y, 1.0)
        else:  # forward
            # 正向公式: s = ALT_y / ALT_ref
            s_raw = np.where(valid_mask, alt_y / ALT_REF, 1.0)

    # 处理无效值
    s_raw = np.where(np.isnan(s_raw) | np.isinf(s_raw), 1.0, s_raw)

    # 截断到合理范围
    if clip_low is not None and clip_high is not None:
        s_raw = np.clip(s_raw, float(clip_low), float(clip_high))

    # 放大：delta * alpha
    permafrost_mask_y = valid_mask
    return s_raw, alt_y, permafrost_mask_y


def _get_alpha_for_param(alpha_cfg, param):
    if isinstance(alpha_cfg, dict):
        if param in alpha_cfg:
            return float(alpha_cfg[param])
        if 'default' in alpha_cfg:
            return float(alpha_cfg['default'])
        return 1.0
    return float(alpha_cfg)


def apply_scheme_params(scheme, year, base_shape, return_stats=False):
    """根据方案计算修正后的参数栅格"""
    fc_grid = np.full(base_shape, V6_PARAMS['FC'], dtype=np.float64)
    perc_grid = np.full(base_shape, V6_PARAMS['PERC'], dtype=np.float64)
    k2_grid = np.full(base_shape, V6_PARAMS['K2'], dtype=np.float64)
    k_grid = np.full(base_shape, V6_PARAMS['K'], dtype=np.float64)

    use_alt = bool(scheme.get('use_alt', False))

    # Optional: activate coupling only for a subset of years (e.g., post-2017 regime shift)
    start_year = scheme.get('start_year', None)
    end_year = scheme.get('end_year', None)
    apply_coupling = use_alt
    if start_year is not None and year < int(start_year):
        apply_coupling = False
    if end_year is not None and year > int(end_year):
        apply_coupling = False

    s_raw, alt_y, permafrost_mask = (None, None, None)
    if use_alt:
        s_raw, alt_y, permafrost_mask = compute_scale_factor(year, scheme)

    s_nl = None
    alpha_cfg = scheme.get('alpha', 1.0)
    power = float(scheme.get('power', 1.0))
    inverse_params = scheme.get('inverse_params', [])

    if apply_coupling and (s_raw is not None) and (permafrost_mask is not None):
        with np.errstate(all='ignore'):
            s_nl = np.power(s_raw, power) if abs(power - 1.0) > 1e-12 else s_raw

        # Optional: year-ramp for coupling strength (0..1), applied to alpha
        year_weight = 1.0
        ramp = scheme.get('year_ramp', None)
        if ramp is not None:
            y0, y1 = int(ramp[0]), int(ramp[1])
            if y1 == y0:
                year_weight = 1.0 if year >= y1 else 0.0
            else:
                year_weight = (year - y0) / float(y1 - y0)
                year_weight = float(np.clip(year_weight, 0.0, 1.0))

        for param in scheme.get('params', []):
            alpha_p = _get_alpha_for_param(alpha_cfg, param) * year_weight

            # Optional: apply different coupling strength by elevation zones
            # zone_weights example: {'high': 1.5, 'mid': 1.0, 'low': 0.7}
            zone_weights = scheme.get('zone_weights', None)
            if zone_weights is not None and ZONE_MID is not None and ZONE_HIGH is not None:
                w_high = float(zone_weights.get('high', 1.0))
                w_mid = float(zone_weights.get('mid', 1.0))
                w_low = float(zone_weights.get('low', 1.0))
                w_grid = np.ones_like(s_nl, dtype=np.float64) * w_low
                w_grid[ZONE_MID] = w_mid
                w_grid[ZONE_HIGH] = w_high
                alpha_eff = alpha_p * w_grid
            else:
                alpha_eff = alpha_p

            amp_mode = scheme.get('amp_mode', 'linear')
            if amp_mode == 'exp':
                s_factor = np.power(s_nl, alpha_eff)
            else:
                s_factor = 1.0 + (s_nl - 1.0) * alpha_eff

            factor_clip = scheme.get('factor_clip', None)
            if factor_clip is not None:
                lo, hi = float(factor_clip[0]), float(factor_clip[1])
                s_factor = np.clip(s_factor, lo, hi)

            if param == 'FC':
                if param in inverse_params:
                    fc_grid[permafrost_mask] = np.clip(
                        V6_PARAMS['FC'] / s_factor[permafrost_mask],
                        *PARAM_BOUNDS['FC']
                    )
                else:
                    fc_grid[permafrost_mask] = np.clip(
                        V6_PARAMS['FC'] * s_factor[permafrost_mask],
                        *PARAM_BOUNDS['FC']
                    )

            elif param == 'PERC':
                if param in inverse_params:
                    perc_grid[permafrost_mask] = np.clip(
                        V6_PARAMS['PERC'] / s_factor[permafrost_mask],
                        *PARAM_BOUNDS['PERC']
                    )
                else:
                    perc_grid[permafrost_mask] = np.clip(
                        V6_PARAMS['PERC'] * s_factor[permafrost_mask],
                        *PARAM_BOUNDS['PERC']
                    )

            elif param == 'K2':
                if param in inverse_params:
                    k2_grid[permafrost_mask] = np.clip(
                        V6_PARAMS['K2'] / s_factor[permafrost_mask],
                        *PARAM_BOUNDS['K2']
                    )
                else:
                    k2_grid[permafrost_mask] = np.clip(
                        V6_PARAMS['K2'] * s_factor[permafrost_mask],
                        *PARAM_BOUNDS['K2']
                    )

            elif param in ['K0', 'K']:
                if param in inverse_params:
                    k_grid[permafrost_mask] = np.clip(
                        V6_PARAMS['K'] / s_factor[permafrost_mask],
                        *PARAM_BOUNDS['K']
                    )
                else:
                    k_grid[permafrost_mask] = np.clip(
                        V6_PARAMS['K'] * s_factor[permafrost_mask],
                        *PARAM_BOUNDS['K']
                    )

    # If coupling is disabled for this year, still compute s_nl for interpretability (no effect on params)
    if (s_nl is None) and (s_raw is not None):
        with np.errstate(all='ignore'):
            s_nl = np.power(s_raw, power) if abs(power - 1.0) > 1e-12 else s_raw

    stats = None
    if return_stats:
        stats = {'year': int(year)}

        if (alt_y is not None) and (permafrost_mask is not None) and np.any(permafrost_mask):
            alt_vals = alt_y[permafrost_mask]
            stats['alt_mean_m'] = float(np.nanmean(alt_vals))
            stats['alt_p05_m'] = float(np.nanpercentile(alt_vals, 5))
            stats['alt_p95_m'] = float(np.nanpercentile(alt_vals, 95))

            s_vals = s_raw[permafrost_mask] if s_raw is not None else None
            if s_vals is not None:
                stats['s_raw_mean'] = float(np.nanmean(s_vals))
                stats['s_raw_p05'] = float(np.nanpercentile(s_vals, 5))
                stats['s_raw_p95'] = float(np.nanpercentile(s_vals, 95))

            snl_vals = s_nl[permafrost_mask] if s_nl is not None else None
            if snl_vals is not None:
                stats['s_nl_mean'] = float(np.nanmean(snl_vals))

        # Parameter means (report over冻土掩膜；若无掩膜则全域均值)
        if (permafrost_mask is not None) and np.any(permafrost_mask):
            stats['fc_mean'] = float(np.nanmean(fc_grid[permafrost_mask]))
            stats['perc_mean'] = float(np.nanmean(perc_grid[permafrost_mask]))
            stats['k2_mean'] = float(np.nanmean(k2_grid[permafrost_mask]))
            stats['k_mean'] = float(np.nanmean(k_grid[permafrost_mask]))
        else:
            stats['fc_mean'] = float(np.nanmean(fc_grid))
            stats['perc_mean'] = float(np.nanmean(perc_grid))
            stats['k2_mean'] = float(np.nanmean(k2_grid))
            stats['k_mean'] = float(np.nanmean(k_grid))

        stats['coupling_applied'] = bool(apply_coupling and (s_raw is not None) and (permafrost_mask is not None))

    return (fc_grid, perc_grid, k2_grid, k_grid, stats) if return_stats else (fc_grid, perc_grid, k2_grid, k_grid)


# =============================================================================
# 模拟运行
# =============================================================================
def run_scheme(scheme_name, scheme):
    """运行单个方案"""
    print(f"\n--- {scheme_name}: {scheme['description']} ---")

    TT = V6_PARAMS['TT']
    RFCF, SFCF = V6_PARAMS['RFCF'], V6_PARAMS['SFCF']
    CWH, CFR = V6_PARAMS['CWH'], V6_PARAMS['CFR']
    FC, BETA, LP = V6_PARAMS['FC'], V6_PARAMS['BETA'], V6_PARAMS['LP']
    K, K1, K2 = V6_PARAMS['K'], V6_PARAMS['K1'], V6_PARAMS['K2']
    PERC = V6_PARAMS['PERC']
    K_MUSK = V6_PARAMS['K_MUSK']
    ICE_FACTOR = V6_PARAMS['ICE_FACTOR']
    CFMAX_mid, CFMAX_high = V6_PARAMS['CFMAX_mid'], V6_PARAMS['CFMAX_high']

    cfmax_grid = np.zeros((PREC_3D.shape[0], PREC_3D.shape[1]), dtype=np.float64)
    cfmax_grid[ZONE_MID] = CFMAX_mid
    cfmax_grid[ZONE_HIGH] = CFMAX_high
    cfmax_grid[~(ZONE_MID | ZONE_HIGH)] = CFMAX_mid

    par_base = np.array([
        TT, RFCF, SFCF, 0.0, CWH, CFR,
        FC, BETA, FIXED['ALPHA'], LP, K, K1, K2, FIXED['UZL'], PERC
    ], dtype=np.float64)

    n_cells = len(VALID_CELLS)
    init_st_grid = np.zeros((n_cells, 5), dtype=np.float64)
    for i in range(n_cells):
        init_st_grid[i, :] = INIT_ST

    all_dates = pd.date_range(WARMUP_START, SIM_END, freq='D')
    years = sorted(set(all_dates.year))

    q_total_all, q_rain_all, q_snow_all, q_ice_all = [], [], [], []
    yearly_stats = []

    for year in years:
        year_mask = all_dates.year == year
        year_indices = np.where(year_mask)[0]
        n_days = len(year_indices)
        start_idx, end_idx = year_indices[0], year_indices[-1] + 1

        prec_year = PREC_3D[:, :, start_idx:end_idx].copy()
        temp_year = TEMP_3D[:, :, start_idx:end_idx].copy()
        et_year = ET_3D[:, :, start_idx:end_idx].copy()
        ll_temp_year = LL_TEMP_3D[:, :, start_idx:end_idx].copy()

        fc_grid, perc_grid, k2_grid, k_grid, stat = apply_scheme_params(
            scheme, year, (PREC_3D.shape[0], PREC_3D.shape[1]), return_stats=True
        )
        if stat is not None:
            yearly_stats.append(stat)

        q_total, q_rain, q_snow, q_ice, end_st_grid = run_all_cells_scheme(
            prec_year, temp_year, et_year, ll_temp_year, par_base,
            cfmax_grid, fc_grid, perc_grid, k2_grid, k_grid, init_st_grid,
            VALID_CELLS, CELL_SCALE, GLACIER_MASK, ICE_FACTOR
        )

        init_st_grid = end_st_grid.copy()
        q_total_all.append(q_total)
        q_rain_all.append(q_rain)
        q_snow_all.append(q_snow)
        q_ice_all.append(q_ice)

    q_total_full = np.concatenate(q_total_all)
    q_rain_full = np.concatenate(q_rain_all)
    q_snow_full = np.concatenate(q_snow_all)
    q_ice_full = np.concatenate(q_ice_all)

    q_sim = muskingum_route(q_total_full, K_MUSK, FIXED['X_MUSK'])
    q_rain_routed = muskingum_route(q_rain_full, K_MUSK, FIXED['X_MUSK'])
    q_snow_routed = muskingum_route(q_snow_full, K_MUSK, FIXED['X_MUSK'])
    q_ice_routed = muskingum_route(q_ice_full, K_MUSK, FIXED['X_MUSK'])

    q_sim = q_sim[WARMUP_DAYS:]
    q_rain_routed = q_rain_routed[WARMUP_DAYS:]
    q_snow_routed = q_snow_routed[WARMUP_DAYS:]
    q_ice_routed = q_ice_routed[WARMUP_DAYS:]

    return q_sim, q_rain_routed, q_snow_routed, q_ice_routed, yearly_stats


def calculate_metrics(q_sim, q_obs):
    """计算性能指标"""
    calib_mask = (SIM_DATES >= pd.to_datetime(CALIB_START)) & (SIM_DATES <= pd.to_datetime(CALIB_END))
    valid_mask = (SIM_DATES >= pd.to_datetime(VALID_START)) & (SIM_DATES <= pd.to_datetime(VALID_END))

    q_sim_cal, q_obs_cal = q_sim[calib_mask], q_obs[calib_mask]
    q_sim_val, q_obs_val = q_sim[valid_mask], q_obs[valid_mask]

    nse_cal = float(nse_calc(q_obs_cal, q_sim_cal))
    nse_val = float(nse_calc(q_obs_val, q_sim_val))

    den_cal = np.nansum(q_obs_cal)
    pbias_cal = 0.0 if den_cal <= 0 else 100.0 * np.nansum(q_sim_cal - q_obs_cal) / den_cal

    den_val = np.nansum(q_obs_val)
    pbias_val = 0.0 if den_val <= 0 else 100.0 * np.nansum(q_sim_val - q_obs_val) / den_val

    valid_cal = ~np.isnan(q_obs_cal) & ~np.isnan(q_sim_cal)
    rmse_cal = float(np.sqrt(np.mean((q_obs_cal[valid_cal] - q_sim_cal[valid_cal])**2)))

    valid_val = ~np.isnan(q_obs_val) & ~np.isnan(q_sim_val)
    rmse_val = float(np.sqrt(np.mean((q_obs_val[valid_val] - q_sim_val[valid_val])**2)))

    return {
        'nse_cal': nse_cal, 'nse_val': nse_val,
        'pbias_cal': pbias_cal, 'pbias_val': pbias_val,
        'rmse_cal': rmse_cal, 'rmse_val': rmse_val,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Experiment 2: multi-scheme ALT coupling comparison")
    parser.add_argument(
        "--only",
        type=str,
        default="",
        help="Comma-separated scheme names to run (default: run all). Example: --only baseline,S10_fc_perc_k2_a5",
    )
    parser.add_argument(
        "--alt-ref",
        type=str,
        default="",
        help="ALT reference raster filename (in data/alt/scale_factors) or an absolute path. "
             "Default: use ALT_ref_2006_2017.tif if it exists, else ALT_ref.tif.",
    )
    args = parser.parse_args()

    # Choose ALT_ref file (avoid leakage by preferring calibration-period reference when available)
    global ALT_REF_FILE
    if args.alt_ref.strip():
        ALT_REF_FILE = args.alt_ref.strip()
    else:
        cand = os.path.join(SF_DIR, "ALT_ref_2006_2017.tif")
        ALT_REF_FILE = "ALT_ref_2006_2017.tif" if os.path.exists(cand) else "ALT_ref.tif"

    print("\n" + "="*70)
    print("多年冻土-HBV耦合模型 多方案对比实验")
    print("="*70)

    start_time = time.time()

    # 加载数据
    load_all_data()

    # JIT预热
    print("\nNumba JIT预热...")
    dummy_prec = np.random.rand(5, 5, 10).astype(np.float32)
    dummy_temp = np.random.rand(5, 5, 10).astype(np.float32) * 10
    dummy_et = np.random.rand(5, 5, 10).astype(np.float32)
    dummy_ll = np.random.rand(5, 5, 10).astype(np.float32) * 5
    dummy_par = np.array([0.0, 1.0, 1.0, 3.0, 0.05, 0.05, 200.0, 1.0, 0.0, 0.5, 0.2, 0.05, 0.005, 30.0, 2.0], dtype=np.float64)
    dummy_st = np.array([0.0, 5.0, 0.0, 0.0, 0.0], dtype=np.float64)
    _ = hbv_cell_full(dummy_prec[0,0,:], dummy_temp[0,0,:], dummy_et[0,0,:], dummy_ll[0,0,:], dummy_par, dummy_st, True, 1.5)

    dummy_cells = np.array([[0,0], [1,1]], dtype=np.int64)
    dummy_scale = np.array([1.0, 1.0], dtype=np.float64)
    dummy_glacier = np.zeros((5,5), dtype=np.bool_)
    dummy_cfmax = np.ones((5,5), dtype=np.float64) * 4.0
    dummy_fc = np.ones((5,5), dtype=np.float64) * 200.0
    dummy_perc = np.ones((5,5), dtype=np.float64) * 2.0
    dummy_k2 = np.ones((5,5), dtype=np.float64) * 0.001
    dummy_k = np.ones((5,5), dtype=np.float64) * 0.07
    dummy_init_st = np.zeros((2, 5), dtype=np.float64)
    _ = run_all_cells_scheme(dummy_prec, dummy_temp, dummy_et, dummy_ll, dummy_par,
                             dummy_cfmax, dummy_fc, dummy_perc, dummy_k2, dummy_k, dummy_init_st,
                             dummy_cells, dummy_scale, dummy_glacier, 1.5)
    print("JIT编译完成!")

    # 运行所有方案
    print("\n" + "="*70)
    print("运行多方案对比实验")
    print("="*70)

    schemes_to_run = SCHEMES
    if args.only:
        only_names = [s.strip() for s in args.only.split(",") if s.strip()]
        schemes_to_run = {k: SCHEMES[k] for k in only_names if k in SCHEMES}
        missing = [k for k in only_names if k not in SCHEMES]
        if missing:
            print(f"[WARN] Unknown scheme names skipped: {missing}")

    results = {}
    for scheme_name, scheme in schemes_to_run.items():
        q_sim, q_rain, q_snow, q_ice, yearly_stats = run_scheme(scheme_name, scheme)
        metrics_objective = calculate_metrics(q_sim, Q_OBS_OBJ_FULL)
        metrics_full = calculate_metrics(q_sim, Q_OBS_FULL)

        results[scheme_name] = {
            'description': scheme['description'],
            'metrics': metrics_full,
            'metrics_objective_window': metrics_objective,
            'metrics_full_period': metrics_full,
            'q_sim': q_sim,
            'q_rain': q_rain,
            'q_snow': q_snow,
            'q_ice': q_ice,
            'yearly_stats': yearly_stats,
        }

        print(
            f"   Objective NSE(val)={metrics_objective['nse_val']:.4f}, "
            f"Full NSE(val)={metrics_full['nse_val']:.4f}, "
            f"Full PBIAS(val)={metrics_full['pbias_val']:.1f}%"
        )

    # 结果汇总
    print("\n" + "="*70)
    print("方案对比结果汇总")
    print("="*70)
    print(f"\n{'方案':<20} {'描述':<35} {'NSE(cal)':<10} {'NSE(val)':<10} {'PBIAS(cal)':<12} {'PBIAS(val)':<12}")
    print("-" * 100)

    best_nse_val_obj = -999
    best_scheme_obj = None
    best_nse_val_full = -999
    best_scheme_full = None

    for scheme_name, result in results.items():
        m_obj = result['metrics_objective_window']
        m_full = result['metrics_full_period']
        print(f"{scheme_name:<20} {result['description'][:33]:<35} {m_full['nse_cal']:.4f}     {m_full['nse_val']:.4f}     "
              f"{m_full['pbias_cal']:>+6.1f}%      {m_full['pbias_val']:>+6.1f}%")

        if m_obj['nse_val'] > best_nse_val_obj:
            best_nse_val_obj = m_obj['nse_val']
            best_scheme_obj = scheme_name
        if m_full['nse_val'] > best_nse_val_full:
            best_nse_val_full = m_full['nse_val']
            best_scheme_full = scheme_name

    print("-" * 100)
    print(f"\n最佳方案（目标窗口验证期NSE最高）: {best_scheme_obj}")
    print(f"   {results[best_scheme_obj]['description']}")
    print(f"   验证期 NSE = {best_nse_val_obj:.4f}")
    print(f"\n最佳方案（全时段验证期NSE最高）: {best_scheme_full}")
    print(f"   {results[best_scheme_full]['description']}")
    print(f"   验证期 NSE = {best_nse_val_full:.4f}")

    # 保存结果
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(RESULTS_DIR, f"exp2_multi_scheme_{run_id}")
    os.makedirs(run_dir, exist_ok=True)

    # 保存汇总
    summary = {
        'run_id': run_id,
        'run_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'alt_ref_file': ALT_REF_FILE,
        'obs_mode': {
            'override': OBS_MODE_OVERRIDE,
            'detected': OBS_EVAL_DESC,
            'objective_window': 'may_oct' if OBS_EVAL_DESC.startswith('mayoct') else 'full_year',
            'comparison_metrics': 'full_period',
        },
        'best_scheme': best_scheme_full,
        'best_scheme_objective_window': best_scheme_obj,
        'best_scheme_full_period': best_scheme_full,
        'schemes': {}
    }

    for scheme_name, result in results.items():
        summary['schemes'][scheme_name] = {
            'description': result['description'],
            'metrics': result['metrics'],
            'metrics_objective_window': result['metrics_objective_window'],
            'metrics_full_period': result['metrics_full_period'],
        }

        # 保存时序
        df = pd.DataFrame({
            'date': SIM_DATES,
            'q_sim': result['q_sim'],
            'q_rain': result['q_rain'],
            'q_snow': result['q_snow'],
            'q_ice': result['q_ice'],
            'q_obs': Q_OBS_FULL,
            'q_obs_objective': Q_OBS_OBJ_FULL,
        })
        df.to_csv(os.path.join(run_dir, f'{scheme_name}_timeseries.csv'), index=False)

        # 保存逐年耦合诊断（ALT/尺度因子/参数均值）
        if result.get('yearly_stats'):
            pd.DataFrame(result['yearly_stats']).to_csv(
                os.path.join(run_dir, f'{scheme_name}_param_stats.csv'),
                index=False
            )

    with open(os.path.join(run_dir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    elapsed = time.time() - start_time
    print(f"\n结果已保存到: {run_dir}")
    print(f"运行完成! 耗时: {elapsed/60:.1f} 分钟")

    return results


if __name__ == "__main__":
    main()


