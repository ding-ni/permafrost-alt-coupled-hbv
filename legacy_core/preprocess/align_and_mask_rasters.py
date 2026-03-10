# -*- coding: utf-8 -*-
"""
步骤3：数据对齐与流域掩膜

将所有气象栅格对齐到 DEM 的坐标系、范围和分辨率，
然后裁剪到流域边界（masked）

输入：
- data/gis/dem_1km.tif                    参考DEM
- data/gis/basin.shp             流域边界
- data/raw/precipitation/daily/*.tif      原始降水
- data/raw/temperature/daily/*.tif        原始温度
- data/raw/evaporation/daily/*.tif        原始蒸散发

输出：
- data/aligned_masked/prec/*.tif         对齐并裁剪的降水
- data/aligned_masked/temp/*.tif         对齐并裁剪的温度
- data/aligned_masked/evap/*.tif         对齐并裁剪的蒸散发

修复日志：
- 2026-01-09: 修复文件名按字典序排序的bug，改为按日期排序
- 2026-01-09: 整合流域掩膜步骤，直接输出到aligned_masked
"""

import sys
import os
import glob
import re
from datetime import datetime, timedelta

sys.path.insert(0, r".\Hapi-main\Hapi-main\src")

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.enums import Resampling as ResamplingEnum
from rasterio.mask import mask
import geopandas as gpd

# ============================================================
# 路径配置
# ============================================================
PROJECT_ROOT = r"."

DEM_FILE = os.path.join(PROJECT_ROOT, "data", "gis", "dem_1km.tif")
BASIN_SHP = os.path.join(PROJECT_ROOT, "data", "gis", "basin.shp")

INPUT_DIRS = {
    "prec": os.path.join(PROJECT_ROOT, "data", "raw", "precipitation", "daily"),
    "temp": os.path.join(PROJECT_ROOT, "data", "raw", "temperature", "daily"),
    "evap": os.path.join(PROJECT_ROOT, "data", "raw", "evaporation", "daily"),
}

# 直接输出到 aligned_masked
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "aligned_masked")

# 时间范围配置
START_DATE = datetime(2006, 1, 1)
END_DATE = datetime(2020, 12, 31)
EXPECTED_DAYS = (END_DATE - START_DATE).days + 1  # 5479天

# 文件前缀映射
FILE_PREFIX = {
    "prec": "PREC",
    "temp": "TEMP",
    "evap": "EVAP",
}


# ============================================================
# 对齐函数
# ============================================================
def align_raster_to_dem(input_file, dem_profile, dem_shape):
    """
    将单个栅格对齐到 DEM

    参数:
        input_file: 输入栅格路径
        dem_profile: DEM 的 rasterio profile
        dem_shape: DEM 的形状 (height, width)

    返回:
        对齐后的数据数组
    """
    with rasterio.open(input_file) as src:
        # 准备输出数组
        out_data = np.empty(dem_shape, dtype=np.float32)

        # 重投影
        reproject(
            source=rasterio.band(src, 1),
            destination=out_data,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dem_profile['transform'],
            dst_crs=dem_profile['crs'],
            resampling=Resampling.bilinear
        )

    return out_data


def apply_basin_mask(data, basin_shape, dem_profile):
    """
    应用流域边界掩膜

    参数:
        data: 输入数据数组
        basin_shape: 流域边界shape
        dem_profile: DEM profile

    返回:
        掩膜后的数据
    """
    # 创建临时栅格用于掩膜
    temp_tif = "temp_mask.tif"

    # 保存临时文件
    temp_profile = dem_profile.copy()
    temp_profile.update(dtype=np.float32, nodata=-9999)

    with rasterio.open(temp_tif, 'w', **temp_profile) as dst:
        dst.write(data.astype(np.float32), 1)

    # 应用掩膜
    with rasterio.open(temp_tif) as src:
        out_image, out_transform = mask(src, basin_shape, crop=False)
        out_data = out_image[0]

    # 删除临时文件
    os.remove(temp_tif)

    # 将掩膜外的值设为nodata
    out_data[out_data == dem_profile['nodata']] = -9999

    return out_data


# ============================================================
# 处理单个数据类型
# ============================================================
def process_data_type(data_type, dem_profile, dem_shape, basin_shape=None):
    """
    处理单个数据类型（prec/temp/evap）

    参数:
        data_type: 数据类型 ("prec", "temp", "evap")
        dem_profile: DEM profile
        dem_shape: DEM shape
        basin_shape: 流域边界（如果提供，则应用掩膜）

    返回:
        成功处理的文件数
    """
    input_dir = INPUT_DIRS[data_type]
    output_dir = os.path.join(OUTPUT_DIR, data_type)

    if not os.path.exists(input_dir):
        print(f"   [WARN] 输入目录不存在: {input_dir}")
        return 0

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 获取所有 tif 文件
    input_files = glob.glob(os.path.join(input_dir, "*.tif"))

    if not input_files:
        print(f"   [WARN] 未找到 TIF 文件: {input_dir}")
        return 0

    print(f"   找到 {len(input_files)} 个文件...")

    # 收集文件并按日期排序
    file_date_list = []

    for f in input_files:
        fname = os.path.basename(f)
        # 提取日期部分
        match = re.search(r'(\d{4}\.\d{2}\.\d{2})', fname)
        if match:
            date_str = match.group(1)
            try:
                date_obj = datetime.strptime(date_str, "%Y.%m.%d")
                # 检查日期是否在范围内
                if START_DATE <= date_obj <= END_DATE:
                    file_date_list.append((date_obj, f))
            except:
                pass

    # 按日期排序
    file_date_list.sort(key=lambda x: x[0])

    print(f"   有效文件: {len(file_date_list)} 个 (范围: {START_DATE.strftime('%Y.%m.%d')} ~ {END_DATE.strftime('%Y.%m.%d')})")

    # 处理并保存
    count = 0
    prefix = FILE_PREFIX[data_type]

    for i, (date_obj, input_file) in enumerate(file_date_list):
        date_str = date_obj.strftime("%Y.%m.%d")
        output_file = os.path.join(output_dir, f"{i}_{prefix}_{date_str}.tif")

        # 检查是否已存在
        if os.path.exists(output_file):
            count += 1
            continue

        try:
            # 对齐到DEM
            data = align_raster_to_dem(input_file, dem_profile, dem_shape)

            # 应用流域掩膜（如果提供）
            if basin_shape is not None:
                data = apply_basin_mask(data, basin_shape, dem_profile)

            # 保存
            out_profile = dem_profile.copy()
            out_profile.update(dtype=np.float32, count=1, nodata=-9999)

            with rasterio.open(output_file, 'w', **out_profile) as dst:
                dst.write(data.astype(np.float32), 1)

            count += 1

        except Exception as e:
            print(f"   [ERROR] 处理失败 {date_str}: {e}")

        # 进度显示
        if (i + 1) % 500 == 0:
            print(f"      已处理 {i+1}/{len(file_date_list)}")

    print(f"   [OK] 完成: {count}/{len(file_date_list)} 文件")

    # 验证日期连续性
    if len(file_date_list) > 0:
        first_date = file_date_list[0][0]
        last_date = file_date_list[-1][0]
        actual_days = (last_date - first_date).days + 1

        print(f"   [验证] 日期范围: {first_date.strftime('%Y.%m.%d')} ~ {last_date.strftime('%Y.%m.%d')}")
        print(f"   [验证] 实际天数: {len(file_date_list)}, 预期天数: {EXPECTED_DAYS}")

        if len(file_date_list) < EXPECTED_DAYS:
            missing = EXPECTED_DAYS - len(file_date_list)
            print(f"   [WARN] 缺失 {missing} 天的数据!")

    return count


# ============================================================
# 主函数
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("步骤3：数据对齐与流域掩膜")
    print("=" * 60)

    # 检查 DEM
    if not os.path.exists(DEM_FILE):
        print(f"\n[ERROR] DEM 文件不存在: {DEM_FILE}")
        print("   请先将 DEM 放入 data/gis/ 目录")
        sys.exit(1)

    # 读取 DEM 参数
    print(f"\n读取参考 DEM: {DEM_FILE}")
    with rasterio.open(DEM_FILE) as dem:
        dem_profile = dem.profile
        dem_shape = (dem.height, dem.width)

        print(f"   尺寸: {dem.width} x {dem.height}")
        print(f"   CRS: {dem.crs}")
        print(f"   分辨率: {dem.res}")
        print(f"   范围: {dem.bounds}")

    # 读取流域边界（用于掩膜）
    basin_shape = None
    if os.path.exists(BASIN_SHP):
        print(f"\n读取流域边界: {BASIN_SHP}")
        gdf = gpd.read_file(BASIN_SHP)
        basin_shape = [gdf.geometry.iloc[0]]
        print(f"   流域边界已加载，将应用掩膜")
    else:
        print(f"\n[WARN] 流域边界不存在: {BASIN_SHP}")
        print(f"   将不应用流域掩膜")

    # 处理各类数据
    data_types = [
        ("prec", "降水"),
        ("temp", "温度"),
        ("evap", "蒸散发"),
    ]

    for dtype, name in data_types:
        print(f"\n[处理 {name} 数据]")
        process_data_type(dtype, dem_profile, dem_shape, basin_shape)

    print("\n" + "=" * 60)
    print("[OK] 数据对齐与流域掩膜完成！")
    print("=" * 60)
    print(f"\n输出目录: {OUTPUT_DIR}")

    # 统计输出
    print("\n输出文件统计:")
    for dtype, name in data_types:
        output_dir = os.path.join(OUTPUT_DIR, dtype)
        if os.path.exists(output_dir):
            files = glob.glob(os.path.join(output_dir, "*.tif"))
            print(f"   {name}: {len(files)} 个文件")

            # 检查文件命名是否正确
            if len(files) > 0:
                first_file = os.path.basename(sorted(files)[0])
                last_file = os.path.basename(sorted(files)[-1])
                print(f"      起始: {first_file}")
                print(f"      结束: {last_file}")

