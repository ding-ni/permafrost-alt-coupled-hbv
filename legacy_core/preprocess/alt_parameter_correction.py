# -*- coding: utf-8 -*-
"""
Step 9: 基于 ALT 修正 HBV 参数（实验2核心创新）- V6 增强版

根据活动层厚度（ALT）动态修正 FC、PERC 和 K2 参数

理论基础：
- ALT 越厚 -> 活动层土壤容量越大 -> FC 和 PERC 应增大
- ALT 越厚 -> 深层排水通道越畅通 -> K2（基流系数）应增大
- s_y(i) = ALT_ref(i) / ALT_y(i)

V6增强版改进：
1. 非线性放大：使用 s_y^α (α=2.0) 放大尺度因子效应
2. 增加K2参数：基流系数受冻土影响显著
3. 放宽截断：使用1%-99%分位数，保留更多变化
4. 物理约束：限制参数在合理范围内

修正公式：
- s_y_amp = sign(s_y-1) * |s_y-1|^α + 1   # 以1为中心的非线性放大
- FC_y(i) = FC_calibrated * s_y_amp
- PERC_y(i) = PERC_calibrated * s_y_amp
- K2_y(i) = K2_calibrated * s_y_amp

输入：
- data/gis/dem_1km.tif                     参考栅格
- data/alt/basin/ALT_*.tif              裁剪后的流域ALT数据
- V6率定值                                  FC=1401.86, PERC=0.0616, K2=0.00135

输出：
- data/alt/scale_factors/ALT_ref.tif       多年平均参考场
- data/alt/scale_factors/sy_*.tif          逐年尺度因子（原始）
- data/alt/scale_factors/sy_amp_*.tif      逐年尺度因子（放大后）
- data/parameters/alt_adjusted_v6/fc_*.tif    修正后的FC
- data/parameters/alt_adjusted_v6/perc_*.tif  修正后的PERC
- data/parameters/alt_adjusted_v6/k2_*.tif    修正后的K2
"""

import sys
import os
import argparse

sys.path.insert(0, r".\Hapi-main\Hapi-main\src")

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
import matplotlib.pyplot as plt

# ============================================================
# 路径配置
# ============================================================
PROJECT_ROOT = r"."

DEM_FILE = os.path.join(PROJECT_ROOT, "data", "gis", "dem_1km.tif")
ALT_DIR = os.path.join(PROJECT_ROOT, "data", "alt", "basin")
SF_DIR = os.path.join(PROJECT_ROOT, "data", "alt", "scale_factors")
GLOBAL_PARAM_DIR = os.path.join(PROJECT_ROOT, "data", "parameters", "global")
ALT_PARAM_DIR = os.path.join(PROJECT_ROOT, "data", "parameters", "alt_adjusted_v6")  # V6版本输出目录
FIGURES_DIR = os.path.join(PROJECT_ROOT, "results", "figures")

# V6 率定参数（来自 v6_mswep_inline_20260127_182327）
FC_CALIBRATED = 1401.862945   # mm - 土壤最大田间持水量
PERC_CALIBRATED = 0.061616    # mm/day - 渗漏率
K2_CALIBRATED = 0.001354      # 1/day - 基流系数

# 非线性放大指数
AMPLIFICATION_ALPHA = 2.0     # α=2.0 平方放大效应

# 尺度因子公式（与 run_experiment2_multi_scheme.py 保持一致）
# - inverse:  s = ALT_ref / ALT_y
# - forward:  s = ALT_y / ALT_ref
SCALE_FORMULA = "inverse"

# 可选：先做非线性响应 s^p（p=1.0 不启用），再做线性放大
SCALE_POWER = 1.0

# 对某些参数做反向修正（除以尺度因子），例如：K 代表快速出流系数
INVERSE_PARAMS = set()

# 可选：同时生成 K 栅格（快速出流系数）
GENERATE_K = False
K_CALIBRATED = 0.073742
K_MIN, K_MAX = 0.02, 0.5

# 运行标签：非默认配置时用于避免覆盖已有输出
RUN_TAG = ""

# 参数物理约束范围
FC_MIN, FC_MAX = 200.0, 2500.0       # mm
PERC_MIN, PERC_MAX = 0.01, 0.5       # mm/day
K2_MIN, K2_MAX = 0.0002, 0.01        # 1/day

# 时间范围
YEARS = range(2006, 2021)

# ALT_ref 参考场构建的年份范围（默认为 2006-2020）
DEFAULT_REF_START_YEAR = 2006
DEFAULT_REF_END_YEAR = 2020
REF_START_YEAR = DEFAULT_REF_START_YEAR
REF_END_YEAR = DEFAULT_REF_END_YEAR

# 分位数截断阈值（放宽到1%-99%）
QUANTILE_LOW = 0.01
QUANTILE_HIGH = 0.99


# ============================================================
# 读取并重采样ALT到DEM网格
# ============================================================
def read_alt_to_dem_grid(alt_file, dem_profile, dem_shape):
    """读取ALT文件并重采样到DEM网格"""

    with rasterio.open(alt_file) as src:
        alt_data = src.read(1)
        src_nodata = src.nodata if src.nodata is not None else -9999

        # 将nodata转为nan
        alt_data = np.where(alt_data == src_nodata, np.nan, alt_data)

        # 如果尺寸不匹配，需要重采样
        if alt_data.shape != dem_shape:
            # 创建目标数组
            dst_data = np.zeros(dem_shape, dtype=np.float32)

            # 重采样
            reproject(
                source=np.where(np.isnan(alt_data), -9999, alt_data),
                destination=dst_data,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dem_profile['transform'],
                dst_crs=dem_profile['crs'],
                resampling=Resampling.bilinear,
                src_nodata=-9999,
                dst_nodata=-9999
            )

            # 转回nan
            alt_data = np.where(dst_data == -9999, np.nan, dst_data)

    return alt_data.astype(np.float32)


# ============================================================
# 步骤1：计算 ALT 多年平均参考场
# ============================================================
def _get_alt_ref_file():
    """根据 REF_START_YEAR/REF_END_YEAR 返回 ALT_ref 输出文件路径"""
    if REF_START_YEAR == DEFAULT_REF_START_YEAR and REF_END_YEAR == DEFAULT_REF_END_YEAR:
        return os.path.join(SF_DIR, "ALT_ref.tif")
    return os.path.join(SF_DIR, f"ALT_ref_{REF_START_YEAR}_{REF_END_YEAR}.tif")


def calculate_alt_reference(dem_profile, dem_shape):
    """计算 ALT_ref = mean(ALT_y), y in [REF_START_YEAR, REF_END_YEAR]"""

    print("=" * 60)
    print("步骤1：计算 ALT 多年平均参考场")
    print("=" * 60)

    # 收集参考期年份的 ALT
    alt_stack = []

    ref_years = range(int(REF_START_YEAR), int(REF_END_YEAR) + 1)
    for year in ref_years:
        alt_file = os.path.join(ALT_DIR, f"ALT_{year}.tif")

        if not os.path.exists(alt_file):
            print(f"   [WARN] 跳过 {year}: 文件不存在")
            continue

        alt_data = read_alt_to_dem_grid(alt_file, dem_profile, dem_shape)
        alt_stack.append(alt_data)

    if not alt_stack:
        print("[ERROR] 未找到任何 ALT 文件")
        return None

    print(f"   参考期: {REF_START_YEAR}–{REF_END_YEAR}，读取 {len(alt_stack)} 年 ALT 数据")

    # 堆叠为3D数组 [时间, 行, 列]
    alt_array = np.stack(alt_stack, axis=0)

    # 计算多年平均（忽略 NaN）
    with np.errstate(all='ignore'):
        alt_ref = np.nanmean(alt_array, axis=0)

    # 保存
    os.makedirs(SF_DIR, exist_ok=True)
    output_file = _get_alt_ref_file()

    meta = dem_profile.copy()
    meta.update(dtype=np.float32, nodata=-9999)

    alt_ref_save = np.where(np.isnan(alt_ref), -9999, alt_ref)
    with rasterio.open(output_file, "w", **meta) as dst:
        dst.write(alt_ref_save.astype(np.float32), 1)

    # 统计
    valid = alt_ref[~np.isnan(alt_ref)]
    print(f"   [OK] 保存: {output_file}")
    print(f"   ALT_ref 统计: min={valid.min():.2f}, max={valid.max():.2f}, mean={valid.mean():.2f} m")
    print(f"   有效像元: {len(valid)} / {alt_ref.size} ({100*len(valid)/alt_ref.size:.1f}%)")

    return alt_ref


# ============================================================
# 步骤2：计算逐年尺度因子
# ============================================================
def calculate_scale_factors(alt_ref, dem_profile, dem_shape):
    """计算 s_y*(i) = clip(ALT_ref/ALT_y, q_0.05, q_0.95)"""

    print("\n" + "=" * 60)
    print("步骤2：计算逐年尺度因子")
    print("=" * 60)

    # 首先收集所有 s_y 值，用于计算分位数
    all_sy_values = []

    for year in YEARS:
        alt_file = os.path.join(ALT_DIR, f"ALT_{year}.tif")

        if not os.path.exists(alt_file):
            continue

        alt_y = read_alt_to_dem_grid(alt_file, dem_profile, dem_shape)

        # s_y = f(ALT_ref, ALT_y)
        with np.errstate(divide='ignore', invalid='ignore'):
            if SCALE_FORMULA == "forward":
                sy = alt_y / alt_ref
            else:
                sy = alt_ref / alt_y

        # 收集有效值
        valid = ~np.isnan(sy) & ~np.isinf(sy)
        all_sy_values.extend(sy[valid].flatten())

    # 计算分位数阈值
    q_low = np.percentile(all_sy_values, QUANTILE_LOW * 100)
    q_high = np.percentile(all_sy_values, QUANTILE_HIGH * 100)

    print(f"   分位数截断: q_{QUANTILE_LOW} = {q_low:.3f}, q_{QUANTILE_HIGH} = {q_high:.3f}")

    # 准备元数据
    meta = dem_profile.copy()
    meta.update(dtype=np.float32, nodata=-9999)

    # 逐年计算并保存截断后的尺度因子
    for year in YEARS:
        alt_file = os.path.join(ALT_DIR, f"ALT_{year}.tif")

        if not os.path.exists(alt_file):
            print(f"   [WARN] 跳过 {year}")
            continue

        alt_y = read_alt_to_dem_grid(alt_file, dem_profile, dem_shape)

        # s_y = f(ALT_ref, ALT_y)
        with np.errstate(divide='ignore', invalid='ignore'):
            if SCALE_FORMULA == "forward":
                sy = alt_y / alt_ref
            else:
                sy = alt_ref / alt_y

        # 截断到分位数范围
        sy_clipped = np.clip(sy, q_low, q_high)

        # NaN 保持为 NaN
        sy_clipped[np.isnan(alt_y) | np.isnan(alt_ref)] = np.nan

        # 保存
        output_file = os.path.join(SF_DIR, f"sy_{year}{RUN_TAG}.tif")
        sy_save = np.where(np.isnan(sy_clipped), -9999, sy_clipped)
        with rasterio.open(output_file, "w", **meta) as dst:
            dst.write(sy_save.astype(np.float32), 1)

        valid = sy_clipped[~np.isnan(sy_clipped)]
        print(f"   [OK] {year}: mean(s_y*)={valid.mean():.3f}")

    return q_low, q_high


# ============================================================
# 步骤3：生成逐年修正参数（V6增强版 - 非线性放大 + K2）
# ============================================================
def amplify_scale_factor(sy, alpha=AMPLIFICATION_ALPHA):
    """
    以1为中心线性放大偏离量

    s_y_amp = 1 + (s_y - 1) * α

    例如: α=2.0
    - s_y=0.90 -> delta=-0.10 -> s_y_amp = 1 + (-0.10)*2 = 0.80  (变化从-10%放大到-20%)
    - s_y=1.10 -> delta=+0.10 -> s_y_amp = 1 + (+0.10)*2 = 1.20  (变化从+10%放大到+20%)
    - s_y=1.00 -> delta=0     -> s_y_amp = 1.00                   (无变化)
    """
    delta = sy - 1.0
    sy_amp = 1.0 + delta * alpha
    return sy_amp


def generate_adjusted_parameters(alt_ref, dem_profile):
    """
    生成 ALT 修正后的 FC、PERC 和 K2 参数（V6增强版）

    增强特性：
    1. 非线性放大：delta * α (α=2.0) 放大尺度因子偏离量
    2. 增加K2参数：基流系数受冻土影响
    3. 物理约束：限制参数在合理范围内
    """

    print("\n" + "=" * 60)
    print("步骤3：生成逐年修正参数（V6增强版）")
    print("=" * 60)

    print(f"   V6率定基础参数:")
    print(f"   FC_calibrated = {FC_CALIBRATED:.2f} mm")
    print(f"   PERC_calibrated = {PERC_CALIBRATED:.6f} mm/day")
    print(f"   K2_calibrated = {K2_CALIBRATED:.6f} 1/day")
    print(f"\n   非线性放大系数 α = {AMPLIFICATION_ALPHA}")

    # 多年冻土区掩膜 (有ALT数据的区域)
    permafrost_mask = ~np.isnan(alt_ref)
    print(f"   多年冻土区像元: {permafrost_mask.sum()}")
    print(f"   非冻土区像元: {(~permafrost_mask).sum()}")

    os.makedirs(ALT_PARAM_DIR, exist_ok=True)

    # 准备输出元数据
    meta = dem_profile.copy()
    meta.update(dtype=np.float32, nodata=-9999)

    # 统计信息
    all_fc = []
    all_perc = []
    all_k2 = []
    all_k = []

    # 逐年生成修正参数
    for year in YEARS:
        sy_file = os.path.join(SF_DIR, f"sy_{year}{RUN_TAG}.tif")

        if not os.path.exists(sy_file):
            print(f"   [WARN] 跳过 {year}: 尺度因子不存在")
            continue

        with rasterio.open(sy_file) as src:
            sy_star = src.read(1)
            sy_nodata = src.nodata if src.nodata is not None else -9999
            sy_star = np.where(sy_star == sy_nodata, np.nan, sy_star)

        # 可选：先做非线性响应 s^p，再做线性放大
        with np.errstate(all='ignore'):
            sy_nl = np.power(sy_star, SCALE_POWER) if abs(SCALE_POWER - 1.0) > 1e-12 else sy_star
        sy_amp = amplify_scale_factor(sy_nl, AMPLIFICATION_ALPHA)

        # 保存放大后的尺度因子
        sy_amp_file = os.path.join(SF_DIR, f"sy_amp_{year}{RUN_TAG}.tif")
        sy_amp_save = np.where(np.isnan(sy_amp), -9999, sy_amp)
        with rasterio.open(sy_amp_file, "w", **meta) as dst:
            dst.write(sy_amp_save.astype(np.float32), 1)

        # 初始化为率定值（非冻土区保持不变）
        fc_adjusted = np.full(alt_ref.shape, FC_CALIBRATED, dtype=np.float32)
        perc_adjusted = np.full(alt_ref.shape, PERC_CALIBRATED, dtype=np.float32)
        k2_adjusted = np.full(alt_ref.shape, K2_CALIBRATED, dtype=np.float32)
        k_adjusted = None
        if GENERATE_K:
            k_adjusted = np.full(alt_ref.shape, K_CALIBRATED, dtype=np.float32)

        # 修正参数（仅在冻土区）
        valid_mask = permafrost_mask & ~np.isnan(sy_amp)

        # FC修正（应用物理约束）
        if "FC" in INVERSE_PARAMS:
            fc_adjusted[valid_mask] = np.clip(
                FC_CALIBRATED / sy_amp[valid_mask],
                FC_MIN, FC_MAX
            )
        else:
            fc_adjusted[valid_mask] = np.clip(
                FC_CALIBRATED * sy_amp[valid_mask],
                FC_MIN, FC_MAX
            )

        # PERC修正（应用物理约束）
        if "PERC" in INVERSE_PARAMS:
            perc_adjusted[valid_mask] = np.clip(
                PERC_CALIBRATED / sy_amp[valid_mask],
                PERC_MIN, PERC_MAX
            )
        else:
            perc_adjusted[valid_mask] = np.clip(
                PERC_CALIBRATED * sy_amp[valid_mask],
                PERC_MIN, PERC_MAX
            )

        # K2修正（应用物理约束）
        if "K2" in INVERSE_PARAMS:
            k2_adjusted[valid_mask] = np.clip(
                K2_CALIBRATED / sy_amp[valid_mask],
                K2_MIN, K2_MAX
            )
        else:
            k2_adjusted[valid_mask] = np.clip(
                K2_CALIBRATED * sy_amp[valid_mask],
                K2_MIN, K2_MAX
            )

        # K 修正（可选：快速出流系数）
        if GENERATE_K and k_adjusted is not None:
            if "K" in INVERSE_PARAMS:
                k_adjusted[valid_mask] = np.clip(
                    K_CALIBRATED / sy_amp[valid_mask],
                    K_MIN, K_MAX
                )
            else:
                k_adjusted[valid_mask] = np.clip(
                    K_CALIBRATED * sy_amp[valid_mask],
                    K_MIN, K_MAX
                )

        # 保存 FC
        fc_output = os.path.join(ALT_PARAM_DIR, f"fc_{year}{RUN_TAG}.tif")
        fc_save = np.where(np.isnan(fc_adjusted), -9999, fc_adjusted)
        with rasterio.open(fc_output, "w", **meta) as dst:
            dst.write(fc_save.astype(np.float32), 1)

        # 保存 PERC
        perc_output = os.path.join(ALT_PARAM_DIR, f"perc_{year}{RUN_TAG}.tif")
        perc_save = np.where(np.isnan(perc_adjusted), -9999, perc_adjusted)
        with rasterio.open(perc_output, "w", **meta) as dst:
            dst.write(perc_save.astype(np.float32), 1)

        # 保存 K2
        k2_output = os.path.join(ALT_PARAM_DIR, f"k2_{year}{RUN_TAG}.tif")
        k2_save = np.where(np.isnan(k2_adjusted), -9999, k2_adjusted)
        with rasterio.open(k2_output, "w", **meta) as dst:
            dst.write(k2_save.astype(np.float32), 1)

        # 保存 K（可选）
        if GENERATE_K and k_adjusted is not None:
            k_output = os.path.join(ALT_PARAM_DIR, f"k_{year}{RUN_TAG}.tif")
            k_save = np.where(np.isnan(k_adjusted), -9999, k_adjusted)
            with rasterio.open(k_output, "w", **meta) as dst:
                dst.write(k_save.astype(np.float32), 1)

        # 统计
        fc_permafrost = fc_adjusted[permafrost_mask]
        perc_permafrost = perc_adjusted[permafrost_mask]
        k2_permafrost = k2_adjusted[permafrost_mask]
        fc_valid = fc_permafrost[~np.isnan(fc_permafrost)]
        perc_valid = perc_permafrost[~np.isnan(perc_permafrost)]
        k2_valid = k2_permafrost[~np.isnan(k2_permafrost)]

        all_fc.extend(fc_valid)
        all_perc.extend(perc_valid)
        all_k2.extend(k2_valid)

        k_valid = None
        if GENERATE_K and k_adjusted is not None:
            k_permafrost = k_adjusted[permafrost_mask]
            k_valid = k_permafrost[~np.isnan(k_permafrost)]
            all_k.extend(k_valid)

        # 计算放大后的尺度因子统计
        sy_amp_valid = sy_amp[valid_mask]

        msg = (
            f"   [OK] {year}: s_y_amp={sy_amp_valid.mean():.3f}±{sy_amp_valid.std():.3f}, "
            f"FC={fc_valid.mean():.0f}±{fc_valid.std():.0f}, "
            f"PERC={perc_valid.mean():.4f}±{perc_valid.std():.4f}, "
            f"K2={k2_valid.mean():.6f}±{k2_valid.std():.6f}"
        )
        if k_valid is not None and len(k_valid) > 0:
            msg += f", K={k_valid.mean():.6f}±{k_valid.std():.6f}"
        print(msg)

    # 总体统计
    print("\n" + "-" * 60)
    print("   参数变化范围（所有年份冻土区）:")
    all_fc = np.array(all_fc)
    all_perc = np.array(all_perc)
    all_k2 = np.array(all_k2)
    all_k = np.array(all_k) if all_k else None
    print(f"   FC:   {all_fc.min():.0f} ~ {all_fc.max():.0f} mm "
          f"(基准{FC_CALIBRATED:.0f}, 变化{100*(all_fc.min()/FC_CALIBRATED-1):.1f}% ~ {100*(all_fc.max()/FC_CALIBRATED-1):.1f}%)")
    print(f"   PERC: {all_perc.min():.4f} ~ {all_perc.max():.4f} mm/d "
          f"(基准{PERC_CALIBRATED:.4f}, 变化{100*(all_perc.min()/PERC_CALIBRATED-1):.1f}% ~ {100*(all_perc.max()/PERC_CALIBRATED-1):.1f}%)")
    print(f"   K2:   {all_k2.min():.6f} ~ {all_k2.max():.6f} 1/d "
          f"(基准{K2_CALIBRATED:.6f}, 变化{100*(all_k2.min()/K2_CALIBRATED-1):.1f}% ~ {100*(all_k2.max()/K2_CALIBRATED-1):.1f}%)")
    if all_k is not None:
        print(f"   K:    {all_k.min():.6f} ~ {all_k.max():.6f} 1/d "
              f"(基准{K_CALIBRATED:.6f}, 变化{100*(all_k.min()/K_CALIBRATED-1):.1f}% ~ {100*(all_k.max()/K_CALIBRATED-1):.1f}%)")

    return True


# ============================================================
# 可视化
# ============================================================
def visualize_results():
    """可视化 ALT 修正结果"""

    print("\n" + "=" * 60)
    print("生成可视化图表")
    print("=" * 60)

    os.makedirs(FIGURES_DIR, exist_ok=True)

    # 图1: ALT_ref
    alt_ref_file = _get_alt_ref_file()
    if os.path.exists(alt_ref_file):
        with rasterio.open(alt_ref_file) as src:
            alt_ref = src.read(1)
            alt_ref = np.where(alt_ref == src.nodata, np.nan, alt_ref)

        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.imshow(alt_ref, cmap='viridis')
        plt.colorbar(im, ax=ax, label='ALT (m)')
        ax.set_title('Multi-year Mean Active Layer Thickness (ALT_ref)')
        plt.savefig(os.path.join(FIGURES_DIR, "ALT_ref.png"), dpi=150, bbox_inches='tight')
        plt.close()
        print("   [OK] ALT_ref.png")

    # 图2: 尺度因子时间序列
    sy_means = []
    for year in YEARS:
        sy_file = os.path.join(SF_DIR, f"sy_{year}{RUN_TAG}.tif")
        if os.path.exists(sy_file):
            with rasterio.open(sy_file) as src:
                sy = src.read(1)
                sy = np.where(sy == src.nodata, np.nan, sy)
            sy_means.append((year, np.nanmean(sy)))

    if sy_means:
        years, means = zip(*sy_means)
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.plot(years, means, 'o-', linewidth=2, markersize=8)
        ax.axhline(y=1.0, color='r', linestyle='--', label='s_y = 1')
        ax.set_xlabel('Year')
        ax.set_ylabel('Mean Scale Factor s_y*')
        ax.set_title('ALT Scale Factor Interannual Variation')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.savefig(os.path.join(FIGURES_DIR, "scale_factor_timeseries.png"), dpi=150, bbox_inches='tight')
        plt.close()
        print("   [OK] scale_factor_timeseries.png")


# ============================================================
# 主函数
# ============================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate ALT-coupled HBV parameter rasters.")
    parser.add_argument("--formula", choices=["inverse", "forward"], default=SCALE_FORMULA,
                        help="Scale factor formula: inverse=ALT_ref/ALT_y, forward=ALT_y/ALT_ref.")
    parser.add_argument("--alpha", type=float, default=AMPLIFICATION_ALPHA,
                        help="Linear amplification coefficient around 1.0.")
    parser.add_argument("--power", type=float, default=SCALE_POWER,
                        help="Optional nonlinear response: use s^power before amplification.")
    parser.add_argument("--inverse-params", type=str, default="",
                        help="Comma-separated params to invert (divide by scale factor), e.g. 'K'.")
    parser.add_argument("--include-k", action="store_true", help="Also generate K rasters (quick response).")
    parser.add_argument("--ref-start-year", type=int, default=REF_START_YEAR,
                        help="Start year for ALT_ref reference field (default: 2006).")
    parser.add_argument("--ref-end-year", type=int, default=REF_END_YEAR,
                        help="End year for ALT_ref reference field (default: 2020).")
    args = parser.parse_args()

    # Override globals for this run
    SCALE_FORMULA = args.formula
    AMPLIFICATION_ALPHA = float(args.alpha)
    SCALE_POWER = float(args.power)
    INVERSE_PARAMS = {p.strip().upper() for p in args.inverse_params.split(",") if p.strip()}
    GENERATE_K = bool(args.include_k)
    REF_START_YEAR = int(args.ref_start_year)
    REF_END_YEAR = int(args.ref_end_year)

    if REF_END_YEAR < REF_START_YEAR:
        raise ValueError(f"Invalid ref year range: {REF_START_YEAR}–{REF_END_YEAR}")

    def _fmt_float(x: float) -> str:
        s = f"{x:g}"
        return s.replace("-", "m").replace(".", "p")

    # Default run keeps legacy filenames (no tag) to stay compatible with existing pipelines.
    RUN_TAG = ""
    if (
        SCALE_FORMULA != "inverse"
        or abs(AMPLIFICATION_ALPHA - 2.0) > 1e-12
        or abs(SCALE_POWER - 1.0) > 1e-12
        or INVERSE_PARAMS
        or GENERATE_K
        or REF_START_YEAR != DEFAULT_REF_START_YEAR
        or REF_END_YEAR != DEFAULT_REF_END_YEAR
    ):
        parts = [SCALE_FORMULA, f"a{_fmt_float(AMPLIFICATION_ALPHA)}", f"p{_fmt_float(SCALE_POWER)}"]
        if INVERSE_PARAMS:
            parts.append("inv" + "_".join(sorted(INVERSE_PARAMS)))
        if GENERATE_K:
            parts.append("K")
        if REF_START_YEAR != DEFAULT_REF_START_YEAR or REF_END_YEAR != DEFAULT_REF_END_YEAR:
            parts.append(f"ref{REF_START_YEAR}_{REF_END_YEAR}")
        RUN_TAG = "_" + "_".join(parts)

    print(f"[CFG] ref_years={REF_START_YEAR}-{REF_END_YEAR}, formula={SCALE_FORMULA}, alpha={AMPLIFICATION_ALPHA:g}, power={SCALE_POWER:g}, "
          f"inverse_params={sorted(INVERSE_PARAMS)}, include_k={GENERATE_K}, tag='{RUN_TAG}'")

    print("Step 9: 目标流域 ALT 参数修正（V6版本）")
    print("=" * 60)
    print(f"基础参数: FC = {FC_CALIBRATED:.2f} mm, PERC = {PERC_CALIBRATED:.6f} mm/d")

    # 读取DEM作为参考
    if not os.path.exists(DEM_FILE):
        print(f"[ERROR] DEM文件不存在: {DEM_FILE}")
        sys.exit(1)

    with rasterio.open(DEM_FILE) as dem:
        dem_data = dem.read(1)
        dem_profile = dem.profile.copy()
        dem_shape = (dem.height, dem.width)

    print(f"参考DEM: {dem_shape}")

    # 检查 ALT 数据
    alt_files = [f for f in os.listdir(ALT_DIR) if f.endswith('.tif')] if os.path.exists(ALT_DIR) else []

    if not alt_files:
        print(f"\n[ERROR] ALT 数据不存在: {ALT_DIR}")
        print("   请先运行 01b_clip_alt_data.py 裁剪 ALT 数据")
        sys.exit(1)

    print(f"找到 {len(alt_files)} 个 ALT 文件")

    # 步骤1: 计算参考场
    alt_ref = calculate_alt_reference(dem_profile, dem_shape)

    if alt_ref is None:
        sys.exit(1)

    # 步骤2: 计算尺度因子
    calculate_scale_factors(alt_ref, dem_profile, dem_shape)

    # 步骤3: 生成修正参数
    generate_adjusted_parameters(alt_ref, dem_profile)

    # 可视化
    visualize_results()

    print("\n" + "=" * 60)
    print("[OK] ALT 参数修正完成！")
    print("=" * 60)
    print(f"\n输出目录:")
    print(f"   尺度因子: {SF_DIR}")
    print(f"   修正参数: {ALT_PARAM_DIR}")

