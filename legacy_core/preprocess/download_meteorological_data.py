# -*- coding: utf-8 -*-
"""
阶段2：下载气象数据

包含：
1. ERA5-Land 温度和蒸散发下载 (需要 cdsapi)
2. MSWEP 降水数据说明

注意：
- ERA5下载需要先注册 CDS 账号并配置 ~/.cdsapirc
- MSWEP需要手动从官网下载
"""

import sys
import os
from datetime import datetime

sys.path.insert(0, r".\Hapi-main\Hapi-main\src")

# ============================================================
# 路径配置
# ============================================================
PROJECT_ROOT = r"."

RAW_TEMP_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "temperature")
RAW_EVAP_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "evaporation")
RAW_PREC_DIR = os.path.join(PROJECT_ROOT, "data", "raw", "precipitation")

# 目标流域大致范围 [North, West, South, East]
# 可根据实际流域边界调整
TUOTUOHE_BBOX = [35.5, 90.5, 33.0, 93.5]

# 时间范围
START_YEAR = 2006
END_YEAR = 2020


# ============================================================
# ERA5-Land 下载配置
# ============================================================
def check_cdsapi_config():
    """检查 CDS API 配置"""
    import os

    # Windows: C:\Users\<username>\.cdsapirc
    # Linux/Mac: ~/.cdsapirc
    home = os.path.expanduser("~")
    config_file = os.path.join(home, ".cdsapirc")

    if os.path.exists(config_file):
        print(f"[OK] CDS API 配置文件存在: {config_file}")
        return True
    else:
        print(f"[ERROR] CDS API 配置文件不存在: {config_file}")
        print("\n请按以下步骤配置：")
        print("1. 访问 https://cds.climate.copernicus.eu/ 注册账号")
        print("2. 登录后访问 https://cds.climate.copernicus.eu/api-how-to")
        print("3. 复制 API key")
        print(f"4. 创建文件 {config_file}，内容如下：")
        print("   url: https://cds.climate.copernicus.eu/api/v2")
        print("   key: <your-uid>:<your-api-key>")
        return False


def extract_if_zip(file_path):
    """如果是ZIP文件则解压"""
    import zipfile
    import shutil

    with open(file_path, 'rb') as f:
        header = f.read(2)

    if header == b'PK':  # ZIP 文件标识
        print(f"      检测到ZIP格式，正在解压...")
        zip_path = file_path + '.zip'
        os.rename(file_path, zip_path)

        with zipfile.ZipFile(zip_path, 'r') as zf:
            names = zf.namelist()
            zf.extractall(os.path.dirname(file_path))

        os.remove(zip_path)

        # 重命名解压出的文件
        for name in names:
            if name.endswith('.nc'):
                extracted = os.path.join(os.path.dirname(file_path), name)
                if os.path.exists(extracted):
                    shutil.move(extracted, file_path)
                    break

        print(f"      解压完成")


def download_era5_temperature(year):
    """下载 ERA5-Land 2m温度数据"""
    import cdsapi

    c = cdsapi.Client()

    output_file = os.path.join(RAW_TEMP_DIR, f"era5_t2m_{year}.nc")

    if os.path.exists(output_file):
        # 检查是否是有效的 NetCDF
        with open(output_file, 'rb') as f:
            header = f.read(2)
        if header != b'PK':  # 不是 ZIP
            print(f"   {year}: 文件已存在，跳过")
            return

    print(f"   {year}: 下载中...")

    c.retrieve(
        'reanalysis-era5-land',
        {
            'variable': '2m_temperature',
            'year': str(year),
            'month': [f'{m:02d}' for m in range(1, 13)],
            'day': [f'{d:02d}' for d in range(1, 32)],
            'time': ['00:00', '06:00', '12:00', '18:00'],  # 每6小时
            'area': TUOTUOHE_BBOX,
            'format': 'netcdf',
        },
        output_file
    )

    # 自动解压 ZIP 文件
    extract_if_zip(output_file)
    print(f"   {year}: [OK]")


def download_era5_evaporation(year):
    """下载 ERA5-Land 蒸散发数据"""
    import cdsapi

    c = cdsapi.Client()

    output_file = os.path.join(RAW_EVAP_DIR, f"era5_evap_{year}.nc")

    if os.path.exists(output_file):
        # 检查是否是有效的 NetCDF
        with open(output_file, 'rb') as f:
            header = f.read(2)
        if header != b'PK':  # 不是 ZIP
            print(f"   {year}: 文件已存在，跳过")
            return

    print(f"   {year}: 下载中...")

    c.retrieve(
        'reanalysis-era5-land',
        {
            'variable': 'total_evaporation',
            'year': str(year),
            'month': [f'{m:02d}' for m in range(1, 13)],
            'day': [f'{d:02d}' for d in range(1, 32)],
            'time': ['00:00', '06:00', '12:00', '18:00'],
            'area': TUOTUOHE_BBOX,
            'format': 'netcdf',
        },
        output_file
    )

    # 自动解压 ZIP 文件
    extract_if_zip(output_file)
    print(f"   {year}: [OK]")


def download_era5_all():
    """下载所有年份的ERA5数据"""
    print("=" * 60)
    print("下载 ERA5-Land 数据")
    print("=" * 60)

    # 确保目录存在
    os.makedirs(RAW_TEMP_DIR, exist_ok=True)
    os.makedirs(RAW_EVAP_DIR, exist_ok=True)

    # 检查配置
    if not check_cdsapi_config():
        return

    try:
        import cdsapi
    except ImportError:
        print("\n[ERROR] cdsapi 未安装")
        print("   pip install cdsapi")
        return

    # 下载温度
    print("\n[1/2] 下载温度数据")
    for year in range(START_YEAR, END_YEAR + 1):
        try:
            download_era5_temperature(year)
        except Exception as e:
            print(f"   [ERROR] {year}: {e}")

    # 下载蒸散发
    print("\n[2/2] 下载蒸散发数据")
    for year in range(START_YEAR, END_YEAR + 1):
        try:
            download_era5_evaporation(year)
        except Exception as e:
            print(f"   [ERROR] {year}: {e}")

    print("\n[OK] ERA5 下载完成！")


# ============================================================
# MSWEP 下载说明
# ============================================================
def print_mswep_instructions():
    """打印 MSWEP 下载说明"""
    print("=" * 60)
    print("MSWEP v2.8 降水数据下载说明")
    print("=" * 60)
    print("""
MSWEP 需要手动下载：

1. 访问 http://www.gloh2o.org/mswep/
2. 注册账号并登录
3. 选择 MSWEP V2.8 产品
4. 下载参数：
   - 时间范围: 2006-01-01 ~ 2020-12-31
   - 时间分辨率: Daily
   - 空间范围: 33°N-35.5°N, 90.5°E-93.5°E (目标流域)
   - 格式: NetCDF

5. 下载后将文件放入:
   {prec_dir}

文件命名格式示例:
   - MSWEP_2006.nc
   - 或 mswep_daily_2006.nc
""".format(prec_dir=RAW_PREC_DIR))


# ============================================================
# 主函数
# ============================================================
if __name__ == "__main__":
    print("气象数据下载工具")
    print("=" * 60)

    # 显示配置
    print(f"\n流域范围 (N, W, S, E): {TUOTUOHE_BBOX}")
    print(f"时间范围: {START_YEAR} - {END_YEAR}")

    # 打印 MSWEP 说明
    print_mswep_instructions()

    # 询问是否下载 ERA5
    print("\n" + "=" * 60)
    response = input("是否开始下载 ERA5-Land 数据? (y/n): ").strip().lower()

    if response == 'y':
        download_era5_all()
    else:
        print("\n跳过 ERA5 下载")
        print("如需下载，请确保已配置 cdsapi 后重新运行")

