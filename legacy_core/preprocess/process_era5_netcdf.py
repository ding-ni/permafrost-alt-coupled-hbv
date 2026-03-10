# -*- coding: utf-8 -*-
"""
阶段2B：ERA5 NetCDF 数据后处理

将 ERA5-Land NetCDF 数据转换为逐日 GeoTIFF 格式

输入：
- data/raw/temperature/era5_t2m_*.nc
- data/raw/evaporation/era5_evap_*.nc

输出：
- data/raw/temperature/daily/*.tif
- data/raw/evaporation/daily/*.tif
"""

import sys
import os

sys.path.insert(0, r".\Hapi-main\Hapi-main\src")

import numpy as np
import xarray as xr
import rasterio
from rasterio.crs import CRS
from datetime import datetime

# ============================================================
# 路径配置
# ============================================================
PROJECT_ROOT = r"."

RAW_TEMP_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "temperature")
RAW_EVAP_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "evaporation")

TEMP_DAILY_DIR = os.path.join(RAW_TEMP_DIR, "daily")
EVAP_DAILY_DIR = os.path.join(RAW_EVAP_DIR, "daily")

YEARS = range(2006, 2021)


# ============================================================
# 处理温度数据
# ============================================================
def process_temperature(year):
    """处理单年温度数据"""

    nc_file = os.path.join(RAW_TEMP_DIR, f"era5_t2m_{year}.nc")

    if not os.path.exists(nc_file):
        print(f"   [WARN] {year}: NetCDF文件不存在")
        return 0

    print(f"   处理 {year} 年温度数据...")

    # 读取 NetCDF
    ds = xr.open_dataset(nc_file)

    # ERA5 温度变量名可能是 't2m' 或 'VAR_2T'
    var_name = 't2m' if 't2m' in ds.data_vars else list(ds.data_vars)[0]
    temp = ds[var_name]

    # 从 Kelvin 转换为 Celsius
    temp_c = temp - 273.15

    # 检查时间维度名称 (可能是 'time' 或 'valid_time')
    time_dim = 'valid_time' if 'valid_time' in temp_c.dims else 'time'

    # 从小时/6小时聚合到日平均
    temp_daily = temp_c.resample({time_dim: '1D'}).mean()

    # 确保输出目录存在
    os.makedirs(TEMP_DAILY_DIR, exist_ok=True)

    # 获取坐标信息
    lons = temp_daily.longitude.values if 'longitude' in temp_daily.coords else temp_daily.lon.values
    lats = temp_daily.latitude.values if 'latitude' in temp_daily.coords else temp_daily.lat.values

    # 计算分辨率
    res_x = abs(lons[1] - lons[0]) if len(lons) > 1 else 0.1
    res_y = abs(lats[1] - lats[0]) if len(lats) > 1 else 0.1

    # 输出计数
    count = 0

    # 获取时间坐标
    time_values = temp_daily[time_dim].values

    # 逐日保存为 GeoTIFF
    for i, time in enumerate(time_values):
        date = np.datetime_as_string(time, unit='D')
        date_str = date.replace('-', '.')

        output_file = os.path.join(TEMP_DAILY_DIR, f"T_{date_str}.tif")

        # 提取该天的数据
        data = temp_daily.isel({time_dim: i}).values

        # 处理 NaN
        data = np.where(np.isnan(data), -9999, data)

        # 创建 transform
        transform = rasterio.transform.from_bounds(
            lons.min() - res_x/2,
            lats.min() - res_y/2,
            lons.max() + res_x/2,
            lats.max() + res_y/2,
            len(lons),
            len(lats)
        )

        # 保存
        with rasterio.open(
            output_file, 'w',
            driver='GTiff',
            height=data.shape[0],
            width=data.shape[1],
            count=1,
            dtype=data.dtype,
            crs=CRS.from_epsg(4326),
            transform=transform,
            nodata=-9999
        ) as dst:
            dst.write(data, 1)

        count += 1

    ds.close()
    print(f"   [OK] {year}: 生成 {count} 个日温度文件")
    return count


# ============================================================
# 处理蒸散发数据
# ============================================================
def process_evaporation(year):
    """处理单年蒸散发数据"""

    nc_file = os.path.join(RAW_EVAP_DIR, f"era5_evap_{year}.nc")

    if not os.path.exists(nc_file):
        print(f"   [WARN] {year}: NetCDF文件不存在")
        return 0

    print(f"   处理 {year} 年蒸散发数据...")

    # 读取 NetCDF
    ds = xr.open_dataset(nc_file)

    # ERA5 蒸散发变量名
    var_name = 'e' if 'e' in ds.data_vars else list(ds.data_vars)[0]
    evap = ds[var_name]

    # ERA5 蒸散发单位是 m，转换为 mm
    # 且 ERA5 蒸散发是负值（向下为正），取绝对值
    evap_mm = np.abs(evap) * 1000

    # 检查时间维度名称 (可能是 'time' 或 'valid_time')
    time_dim = 'valid_time' if 'valid_time' in evap_mm.dims else 'time'

    # 从小时/6小时聚合到日累计
    evap_daily = evap_mm.resample({time_dim: '1D'}).sum()

    # 确保输出目录存在
    os.makedirs(EVAP_DAILY_DIR, exist_ok=True)

    # 获取坐标信息
    lons = evap_daily.longitude.values if 'longitude' in evap_daily.coords else evap_daily.lon.values
    lats = evap_daily.latitude.values if 'latitude' in evap_daily.coords else evap_daily.lat.values

    res_x = abs(lons[1] - lons[0]) if len(lons) > 1 else 0.1
    res_y = abs(lats[1] - lats[0]) if len(lats) > 1 else 0.1

    count = 0

    # 获取时间坐标
    time_values = evap_daily[time_dim].values

    for i, time in enumerate(time_values):
        date = np.datetime_as_string(time, unit='D')
        date_str = date.replace('-', '.')

        output_file = os.path.join(EVAP_DAILY_DIR, f"ET_{date_str}.tif")

        data = evap_daily.isel({time_dim: i}).values
        data = np.where(np.isnan(data), -9999, data)

        transform = rasterio.transform.from_bounds(
            lons.min() - res_x/2,
            lats.min() - res_y/2,
            lons.max() + res_x/2,
            lats.max() + res_y/2,
            len(lons),
            len(lats)
        )

        with rasterio.open(
            output_file, 'w',
            driver='GTiff',
            height=data.shape[0],
            width=data.shape[1],
            count=1,
            dtype=np.float32,
            crs=CRS.from_epsg(4326),
            transform=transform,
            nodata=-9999
        ) as dst:
            dst.write(data.astype(np.float32), 1)

        count += 1

    ds.close()
    print(f"   [OK] {year}: 生成 {count} 个日蒸散发文件")
    return count


# ============================================================
# 主函数
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("ERA5 数据后处理：NetCDF → 逐日 GeoTIFF")
    print("=" * 60)

    # 检查是否有 NetCDF 文件
    temp_files = [f for f in os.listdir(RAW_TEMP_DIR) if f.endswith('.nc')] if os.path.exists(RAW_TEMP_DIR) else []
    evap_files = [f for f in os.listdir(RAW_EVAP_DIR) if f.endswith('.nc')] if os.path.exists(RAW_EVAP_DIR) else []

    print(f"\n找到温度 NetCDF 文件: {len(temp_files)}")
    print(f"找到蒸散发 NetCDF 文件: {len(evap_files)}")

    if not temp_files and not evap_files:
        print("\n[WARN] 未找到 ERA5 NetCDF 文件")
        print("   请先运行 02_download_meteorological_data.py 下载数据")
        sys.exit(0)

    # 处理温度
    if temp_files:
        print("\n[1/2] 处理温度数据")
        total_temp = 0
        for year in YEARS:
            total_temp += process_temperature(year)
        print(f"   共生成 {total_temp} 个温度文件")

    # 处理蒸散发
    if evap_files:
        print("\n[2/2] 处理蒸散发数据")
        total_evap = 0
        for year in YEARS:
            total_evap += process_evaporation(year)
        print(f"   共生成 {total_evap} 个蒸散发文件")

    print("\n[OK] 后处理完成！")

