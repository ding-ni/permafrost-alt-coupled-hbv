# -*- coding: utf-8 -*-
"""
下载 FAO-56 Penman-Monteith 所需的 ERA5 变量

需要的变量:
1. 2m_temperature (已下载)
2. surface_solar_radiation_downwards (太阳辐射)
3. 10m_u_component_of_wind (10m U风速)
4. 10m_v_component_of_wind (10m V风速)
5. 2m_dewpoint_temperature (露点温度)
"""

import os
import sys
import cdsapi
import zipfile
import shutil

# ============================================================
# 配置
# ============================================================
PROJECT_ROOT = r"."
TUOTUOHE_BBOX = [35.5, 90.5, 33.0, 93.5]
START_YEAR = 2006
END_YEAR = 2020

# 输出目录
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
SOLAR_DIR = os.path.join(RAW_DIR, "solar_radiation")
WIND_DIR = os.path.join(RAW_DIR, "wind")
DEWPOINT_DIR = os.path.join(RAW_DIR, "dewpoint")

# 创建目录
for d in [SOLAR_DIR, WIND_DIR, DEWPOINT_DIR]:
    os.makedirs(d, exist_ok=True)


def extract_if_zip(file_path):
    """如果是ZIP文件则解压"""
    with open(file_path, 'rb') as f:
        header = f.read(2)

    if header == b'PK':
        print(f"      解压中...")
        zip_path = file_path + '.zip'
        os.rename(file_path, zip_path)

        with zipfile.ZipFile(zip_path, 'r') as zf:
            names = zf.namelist()
            zf.extractall(os.path.dirname(file_path))

        os.remove(zip_path)

        for name in names:
            if name.endswith('.nc'):
                extracted = os.path.join(os.path.dirname(file_path), name)
                if os.path.exists(extracted):
                    shutil.move(extracted, file_path)
                    break


def download_variable(c, variable, output_dir, prefix, year):
    """下载单个变量单年数据"""
    output_file = os.path.join(output_dir, f"{prefix}_{year}.nc")

    if os.path.exists(output_file):
        size = os.path.getsize(output_file)
        if size > 1000:
            print(f"   {year}: 已存在，跳过")
            return True

    print(f"   {year}: 下载中...")

    try:
        c.retrieve(
            'reanalysis-era5-land',
            {
                'variable': variable,
                'year': str(year),
                'month': [f'{m:02d}' for m in range(1, 13)],
                'day': [f'{d:02d}' for d in range(1, 32)],
                'time': ['00:00', '06:00', '12:00', '18:00'],
                'area': TUOTUOHE_BBOX,
                'format': 'netcdf',
            },
            output_file
        )
        extract_if_zip(output_file)
        print(f"   {year}: [OK]")
        return True
    except Exception as e:
        print(f"   {year}: [ERROR] {e}")
        return False


def main():
    print("=" * 60)
    print("FAO-56 Penman-Monteith 所需变量下载")
    print("=" * 60)

    c = cdsapi.Client()

    # 1. 太阳辐射
    print("\n[1/3] 下载太阳辐射...")
    for year in range(START_YEAR, END_YEAR + 1):
        download_variable(c, 'surface_solar_radiation_downwards',
                         SOLAR_DIR, 'era5_ssrd', year)

    # 2. 风速 (U 和 V 分量)
    print("\n[2/3] 下载风速...")
    for year in range(START_YEAR, END_YEAR + 1):
        # U 分量
        download_variable(c, '10m_u_component_of_wind',
                         WIND_DIR, 'era5_u10', year)
        # V 分量
        download_variable(c, '10m_v_component_of_wind',
                         WIND_DIR, 'era5_v10', year)

    # 3. 露点温度
    print("\n[3/3] 下载露点温度...")
    for year in range(START_YEAR, END_YEAR + 1):
        download_variable(c, '2m_dewpoint_temperature',
                         DEWPOINT_DIR, 'era5_d2m', year)

    print("\n" + "=" * 60)
    print("[OK] FAO-56 变量下载完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()

