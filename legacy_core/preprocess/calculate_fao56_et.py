# -*- coding: utf-8 -*-
"""
计算 FAO-56 Penman-Monteith 参考蒸散发 (ET0)

使用 pyet 库计算标准参考蒸散发

输入：
- data/raw/temperature/era5_t2m_*.nc (2m温度)
- data/raw/solar_radiation/era5_ssrd_*.nc (太阳辐射)
- data/raw/wind/era5_u10_*.nc, era5_v10_*.nc (10m风速)
- data/raw/dewpoint/era5_d2m_*.nc (露点温度)

输出：
- data/raw/evaporation/daily/ET_YYYY.MM.DD.tif (FAO-56 ET0)
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

TEMP_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "temperature")
SOLAR_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "solar_radiation")
WIND_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "wind")
DEWPOINT_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "dewpoint")

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "evaporation", "daily")

YEARS = range(2006, 2021)

# 目标流域平均海拔 (m)
ELEVATION = 4500


# ============================================================
# FAO-56 Penman-Monteith 公式
# ============================================================
def calculate_et0_fao56(t_mean, t_min, t_max, rh_mean, u2, rs, lat, doy, elevation):
    """
    计算 FAO-56 Penman-Monteith 参考蒸散发

    参数:
        t_mean: 日平均温度 (°C)
        t_min: 日最低温度 (°C)
        t_max: 日最高温度 (°C)
        rh_mean: 日平均相对湿度 (%)
        u2: 2m高度风速 (m/s)
        rs: 太阳辐射 (MJ/m²/day)
        lat: 纬度 (度)
        doy: 儒略日 (1-365)
        elevation: 海拔 (m)

    返回:
        et0: 参考蒸散发 (mm/day)
    """
    # 大气压 (kPa)
    P = 101.3 * ((293 - 0.0065 * elevation) / 293) ** 5.26

    # 心理常数 (kPa/°C)
    gamma = 0.665e-3 * P

    # 饱和水汽压 (kPa)
    es_tmax = 0.6108 * np.exp(17.27 * t_max / (t_max + 237.3))
    es_tmin = 0.6108 * np.exp(17.27 * t_min / (t_min + 237.3))
    es = (es_tmax + es_tmin) / 2

    # 实际水汽压 (kPa)
    ea = es * rh_mean / 100

    # 饱和水汽压曲线斜率 (kPa/°C)
    delta = 4098 * (0.6108 * np.exp(17.27 * t_mean / (t_mean + 237.3))) / (t_mean + 237.3) ** 2

    # 日序角 (rad)
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365)

    # 太阳赤纬 (rad)
    delta_s = 0.409 * np.sin(2 * np.pi * doy / 365 - 1.39)

    # 纬度 (rad)
    phi = lat * np.pi / 180

    # 日落时角 (rad)
    omega_s = np.arccos(-np.tan(phi) * np.tan(delta_s))
    omega_s = np.clip(omega_s, 0, np.pi)

    # 外太阳辐射 (MJ/m²/day)
    Gsc = 0.0820  # 太阳常数
    Ra = (24 * 60 / np.pi) * Gsc * dr * (
        omega_s * np.sin(phi) * np.sin(delta_s) +
        np.cos(phi) * np.cos(delta_s) * np.sin(omega_s)
    )
    Ra = np.maximum(Ra, 0)

    # 晴空辐射 (MJ/m²/day)
    Rso = (0.75 + 2e-5 * elevation) * Ra

    # 净短波辐射 (MJ/m²/day)
    albedo = 0.23
    Rns = (1 - albedo) * rs

    # 净长波辐射 (MJ/m²/day)
    sigma = 4.903e-9  # Stefan-Boltzmann 常数
    Rnl = sigma * ((t_max + 273.16) ** 4 + (t_min + 273.16) ** 4) / 2 * \
          (0.34 - 0.14 * np.sqrt(ea)) * \
          (1.35 * np.minimum(rs / np.maximum(Rso, 0.01), 1) - 0.35)

    # 净辐射 (MJ/m²/day)
    Rn = Rns - Rnl

    # 土壤热通量 (日尺度假设为0)
    G = 0

    # FAO-56 Penman-Monteith 公式
    et0 = (0.408 * delta * (Rn - G) + gamma * 900 / (t_mean + 273) * u2 * (es - ea)) / \
          (delta + gamma * (1 + 0.34 * u2))

    # 限制在合理范围
    et0 = np.clip(et0, 0, 15)

    return et0


def calculate_rh_from_dewpoint(t_air, t_dew):
    """从露点温度计算相对湿度"""
    es = 6.112 * np.exp(17.67 * t_air / (t_air + 243.5))
    e = 6.112 * np.exp(17.67 * t_dew / (t_dew + 243.5))
    rh = 100 * e / es
    return np.clip(rh, 0, 100)


def wind_10m_to_2m(u10):
    """将10m风速转换为2m风速"""
    return u10 * 4.87 / np.log(67.8 * 10 - 5.42)


# ============================================================
# 处理单年数据
# ============================================================
def process_year(year):
    """处理单年数据计算ET0"""

    print(f"\n处理 {year} 年...")

    # 检查所有输入文件
    temp_file = os.path.join(TEMP_DIR, f"era5_t2m_{year}.nc")
    solar_file = os.path.join(SOLAR_DIR, f"era5_ssrd_{year}.nc")
    u10_file = os.path.join(WIND_DIR, f"era5_u10_{year}.nc")
    v10_file = os.path.join(WIND_DIR, f"era5_v10_{year}.nc")
    dewpoint_file = os.path.join(DEWPOINT_DIR, f"era5_d2m_{year}.nc")

    missing = []
    for name, f in [("温度", temp_file), ("辐射", solar_file),
                    ("U风速", u10_file), ("V风速", v10_file), ("露点", dewpoint_file)]:
        if not os.path.exists(f):
            missing.append(name)

    if missing:
        print(f"   [WARN] 缺少: {', '.join(missing)}")
        return 0

    # 读取数据
    ds_temp = xr.open_dataset(temp_file)
    ds_solar = xr.open_dataset(solar_file)
    ds_u10 = xr.open_dataset(u10_file)
    ds_v10 = xr.open_dataset(v10_file)
    ds_dew = xr.open_dataset(dewpoint_file)

    # 获取变量
    t2m = ds_temp['t2m'] - 273.15  # K -> °C
    ssrd = ds_solar['ssrd'] / 1e6  # J/m² -> MJ/m²
    u10 = ds_u10['u10']
    v10 = ds_v10['v10']
    d2m = ds_dew['d2m'] - 273.15  # K -> °C

    # 检测时间维度名称 (可能是 'time' 或 'valid_time')
    time_dim = 'valid_time' if 'valid_time' in t2m.dims else 'time'

    # 计算日统计 (使用 dict 式索引以支持动态维度名)
    t_mean = t2m.resample({time_dim: '1D'}).mean()
    t_min = t2m.resample({time_dim: '1D'}).min()
    t_max = t2m.resample({time_dim: '1D'}).max()

    # 处理累积辐射：ERA5-Land每天06 UTC重置，12 UTC达到峰值
    # 正确方法：提取12:00 UTC的值作为日辐射总量
    # 提取每天12:00 UTC的时间步
    times = ssrd[time_dim].values
    hours = [str(t)[11:13] for t in times]
    mask_12utc = [h == '12' for h in hours]

    # 只保留12:00 UTC的数据
    rs_12utc = ssrd.isel({time_dim: mask_12utc})
    # 按日重采样（实际上每天已经只有一个值了）
    rs_daily = rs_12utc.resample({time_dim: '1D'}).mean()

    u10_mean = u10.resample({time_dim: '1D'}).mean()
    v10_mean = v10.resample({time_dim: '1D'}).mean()
    d2m_mean = d2m.resample({time_dim: '1D'}).mean()

    # 计算风速和相对湿度
    wind10 = np.sqrt(u10_mean ** 2 + v10_mean ** 2)
    u2 = wind_10m_to_2m(wind10)  # 转换为2m风速
    rh = calculate_rh_from_dewpoint(t_mean, d2m_mean)

    # 获取坐标
    lons = t_mean.longitude.values if 'longitude' in t_mean.coords else t_mean.lon.values
    lats = t_mean.latitude.values if 'latitude' in t_mean.coords else t_mean.lat.values

    # 创建纬度网格
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    # 确保输出目录存在
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    count = 0
    time_values = t_mean[time_dim].values
    for i, time in enumerate(time_values):
        date = np.datetime_as_string(time, unit='D')
        date_str = date.replace('-', '.')
        doy = int(date[5:7]) * 30 + int(date[8:10])  # 简化的儒略日

        # 提取该天数据
        t_m = t_mean.isel({time_dim: i}).values
        t_mn = t_min.isel({time_dim: i}).values
        t_mx = t_max.isel({time_dim: i}).values
        rh_d = rh.isel({time_dim: i}).values
        u2_d = u2.isel({time_dim: i}).values
        rs_d = rs_daily.isel({time_dim: i}).values

        # 计算 ET0
        et0 = calculate_et0_fao56(
            t_mean=t_m, t_min=t_mn, t_max=t_mx,
            rh_mean=rh_d, u2=u2_d, rs=rs_d,
            lat=lat_grid, doy=doy, elevation=ELEVATION
        )

        # 处理无效值
        et0 = np.where(np.isnan(et0), -9999, et0)

        # 保存
        output_file = os.path.join(OUTPUT_DIR, f"ET_{date_str}.tif")

        res_x = abs(lons[1] - lons[0]) if len(lons) > 1 else 0.1
        res_y = abs(lats[1] - lats[0]) if len(lats) > 1 else 0.1

        transform = rasterio.transform.from_bounds(
            lons.min() - res_x / 2,
            lats.min() - res_y / 2,
            lons.max() + res_x / 2,
            lats.max() + res_y / 2,
            len(lons),
            len(lats)
        )

        with rasterio.open(
            output_file, 'w',
            driver='GTiff',
            height=et0.shape[0],
            width=et0.shape[1],
            count=1,
            dtype=np.float32,
            crs=CRS.from_epsg(4326),
            transform=transform,
            nodata=-9999
        ) as dst:
            dst.write(et0.astype(np.float32), 1)

        count += 1

    # 关闭数据集
    ds_temp.close()
    ds_solar.close()
    ds_u10.close()
    ds_v10.close()
    ds_dew.close()

    print(f"   [OK] 生成 {count} 个 ET0 文件")
    return count


# ============================================================
# 主函数
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("FAO-56 Penman-Monteith 蒸散发计算")
    print("=" * 60)

    # 检查必要文件
    print("\n检查输入数据...")
    for year in YEARS:
        has_all = True
        for d, p in [(TEMP_DIR, "era5_t2m"), (SOLAR_DIR, "era5_ssrd"),
                     (WIND_DIR, "era5_u10"), (WIND_DIR, "era5_v10"),
                     (DEWPOINT_DIR, "era5_d2m")]:
            f = os.path.join(d, f"{p}_{year}.nc")
            if not os.path.exists(f):
                has_all = False
                break
        status = "[OK]" if has_all else "[MISSING]"
        print(f"   {year}: {status}")

    # 处理每年
    print("\n开始计算 FAO-56 ET0...")
    total = 0
    for year in YEARS:
        total += process_year(year)

    print("\n" + "=" * 60)
    print(f"[OK] 完成! 共生成 {total} 个 ET0 文件")
    print(f"输出目录: {OUTPUT_DIR}")
    print("=" * 60)

