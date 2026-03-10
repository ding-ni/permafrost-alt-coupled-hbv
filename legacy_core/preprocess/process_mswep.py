# -*- coding: utf-8 -*-
"""
Step 5-6: MSWEP V2.8 降水数据处理

输入：
- data/raw/precipitation/mswep/*.nc  (从Google Drive下载的MSWEP V2.8 daily文件)
  文件命名格式: YYYYDDD.nc (年+儒略日)

输出：
- data/raw/precipitation/daily/P_YYYY.MM.DD.tif

使用说明：
1. 从Google Drive下载 MSWEP V2.8 -> Past -> Daily 文件夹中的NC文件
2. 只需下载 2006-2020 年的文件 (2006001.nc ~ 2020366.nc)
3. 将文件放入 data/raw/precipitation/mswep/
4. 运行此脚本
"""

import sys
import os
from datetime import datetime, timedelta

sys.path.insert(0, r".\Hapi-main\Hapi-main\src")

import numpy as np
import xarray as xr
import rasterio
from rasterio.crs import CRS

# ============================================================
# 路径配置
# ============================================================
PROJECT_ROOT = r"."

MSWEP_INPUT_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "precipitation", "mswep")
PREC_OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "precipitation", "daily")

# 目标流域范围 [West, South, East, North]
TUOTUOHE_BBOX = {
    'lon_min': 90.5,
    'lon_max': 93.5,
    'lat_min': 33.0,
    'lat_max': 35.5
}

# 时间范围
START_YEAR = 2006
END_YEAR = 2020


# ============================================================
# 儒略日转日期
# ============================================================
def julian_to_date(year, julian_day):
    """将儒略日转换为日期"""
    return datetime(year, 1, 1) + timedelta(days=julian_day - 1)


# ============================================================
# 处理单个MSWEP文件
# ============================================================
def process_mswep_file(nc_file, output_dir):
    """处理单个MSWEP NC文件，提取目标流域区域并保存为GeoTIFF"""

    # 解析文件名获取日期
    basename = os.path.basename(nc_file)
    name = os.path.splitext(basename)[0]

    try:
        year = int(name[:4])
        julian_day = int(name[4:])
        date = julian_to_date(year, julian_day)
        date_str = date.strftime("%Y.%m.%d")
    except:
        print(f"   [WARN] 无法解析文件名: {basename}")
        return False

    output_file = os.path.join(output_dir, f"P_{date_str}.tif")

    # 如果输出文件已存在，跳过
    if os.path.exists(output_file):
        return True

    try:
        # 读取NC文件
        ds = xr.open_dataset(nc_file)

        # MSWEP变量名通常是 'precipitation'
        var_name = 'precipitation' if 'precipitation' in ds.data_vars else list(ds.data_vars)[0]
        prec = ds[var_name]

        # 提取目标流域区域
        # 注意: MSWEP的纬度可能是从大到小排列
        prec_subset = prec.sel(
            lon=slice(TUOTUOHE_BBOX['lon_min'], TUOTUOHE_BBOX['lon_max']),
            lat=slice(TUOTUOHE_BBOX['lat_max'], TUOTUOHE_BBOX['lat_min'])  # 注意顺序
        )

        # 如果有时间维度，取第一个时间步
        if 'time' in prec_subset.dims:
            prec_data = prec_subset.isel(time=0).values
        else:
            prec_data = prec_subset.values

        # 获取坐标
        lons = prec_subset.lon.values
        lats = prec_subset.lat.values

        # 确保纬度是从北到南
        if lats[0] < lats[-1]:
            lats = lats[::-1]
            prec_data = prec_data[::-1, :]

        # 计算分辨率
        res_x = abs(lons[1] - lons[0]) if len(lons) > 1 else 0.1
        res_y = abs(lats[0] - lats[1]) if len(lats) > 1 else 0.1

        # 创建transform
        transform = rasterio.transform.from_bounds(
            lons.min() - res_x/2,
            lats.min() - res_y/2,
            lons.max() + res_x/2,
            lats.max() + res_y/2,
            len(lons),
            len(lats)
        )

        # 处理无效值
        prec_data = np.where(np.isnan(prec_data), -9999, prec_data)
        prec_data = np.where(prec_data < 0, 0, prec_data)  # 降水不能为负

        # 保存为GeoTIFF
        with rasterio.open(
            output_file, 'w',
            driver='GTiff',
            height=prec_data.shape[0],
            width=prec_data.shape[1],
            count=1,
            dtype=np.float32,
            crs=CRS.from_epsg(4326),
            transform=transform,
            nodata=-9999
        ) as dst:
            dst.write(prec_data.astype(np.float32), 1)

        ds.close()
        return True

    except Exception as e:
        print(f"   [ERROR] {basename}: {e}")
        return False


# ============================================================
# 主函数
# ============================================================
def main():
    print("=" * 60)
    print("MSWEP V2.8 降水数据处理")
    print("=" * 60)

    # 检查输入目录
    if not os.path.exists(MSWEP_INPUT_DIR):
        os.makedirs(MSWEP_INPUT_DIR, exist_ok=True)
        print(f"\n[WARN] MSWEP 数据目录不存在，已创建: {MSWEP_INPUT_DIR}")
        print("\n请从 Google Drive 下载 MSWEP V2.8 数据:")
        print("  链接: https://drive.google.com/drive/folders/1Kok05OPVESTpyyan7NafR-2WwuSJ4TO9")
        print("  路径: Past -> Daily")
        print(f"  下载 {START_YEAR}001.nc ~ {END_YEAR}366.nc 文件")
        print(f"  放入: {MSWEP_INPUT_DIR}")
        return

    # 获取所有NC文件
    nc_files = sorted([f for f in os.listdir(MSWEP_INPUT_DIR) if f.endswith('.nc')])

    if not nc_files:
        print(f"\n[WARN] 未找到 NC 文件: {MSWEP_INPUT_DIR}")
        print("\n请从 Google Drive 下载 MSWEP V2.8 数据")
        return

    print(f"\n找到 {len(nc_files)} 个 NC 文件")
    print(f"流域范围: {TUOTUOHE_BBOX}")

    # 创建输出目录
    os.makedirs(PREC_OUTPUT_DIR, exist_ok=True)

    # 处理文件
    success_count = 0
    total = len(nc_files)

    for i, nc_file in enumerate(nc_files):
        nc_path = os.path.join(MSWEP_INPUT_DIR, nc_file)

        if process_mswep_file(nc_path, PREC_OUTPUT_DIR):
            success_count += 1

        # 每100个文件打印进度
        if (i + 1) % 100 == 0 or (i + 1) == total:
            print(f"   进度: {i+1}/{total} ({success_count} 成功)")

    print(f"\n[OK] 处理完成! {success_count}/{total} 文件成功")
    print(f"输出目录: {PREC_OUTPUT_DIR}")

    # 统计输出文件
    output_files = [f for f in os.listdir(PREC_OUTPUT_DIR) if f.endswith('.tif')]
    print(f"生成 {len(output_files)} 个日降水栅格")


if __name__ == "__main__":
    main()

