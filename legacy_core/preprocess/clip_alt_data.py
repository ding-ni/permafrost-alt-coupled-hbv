# -*- coding: utf-8 -*-
"""
阶段1B：裁剪 ALT 数据到目标流域

输入：
- data/gis/basin.shp         流域边界
- data/alt/tibet_plateau/ALT_TP_*.tif  青藏高原ALT数据

输出：
- data/alt/basin/ALT_*.tif  裁剪后的流域ALT数据
"""

import sys
import os

sys.path.insert(0, r".\Hapi-main\Hapi-main\src")

import numpy as np
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import matplotlib.pyplot as plt

# ============================================================
# 路径配置
# ============================================================
PROJECT_ROOT = r"."

BASIN_SHP = os.path.join(PROJECT_ROOT, "data", "gis", "basin.shp")
ALT_INPUT_DIR = os.path.join(PROJECT_ROOT, "data", "alt", "tibet_plateau")
ALT_OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "alt", "basin")

# 时间范围
YEARS = range(2006, 2021)  # 2006-2020


# ============================================================
# 裁剪 ALT 数据
# ============================================================
def clip_alt_to_basin():
    """将青藏高原ALT数据裁剪到流域范围"""

    print("=" * 60)
    print("裁剪 ALT 数据到目标流域")
    print("=" * 60)

    # 检查流域边界
    if not os.path.exists(BASIN_SHP):
        print(f"\n[WARN] 流域边界文件不存在: {BASIN_SHP}")
        print("   请先准备目标流域边界 Shapefile")
        return False

    # 读取流域边界
    print(f"\n读取流域边界: {BASIN_SHP}")
    basin = gpd.read_file(BASIN_SHP)
    print(f"   流域面积: {basin.geometry.area.sum() / 1e6:.2f} km2")
    print(f"   CRS: {basin.crs}")

    # 确保输出目录存在
    os.makedirs(ALT_OUTPUT_DIR, exist_ok=True)

    # 逐年裁剪
    success_count = 0
    for year in YEARS:
        # 尝试不同的文件命名格式
        possible_names = [
            f"ALT_TP_{year}.tif",
            f"ALT_{year}.tif",
            f"alt_{year}.tif",
            f"ALT_TP_{year}_1km.tif",
        ]

        input_file = None
        for name in possible_names:
            path = os.path.join(ALT_INPUT_DIR, name)
            if os.path.exists(path):
                input_file = path
                break

        if input_file is None:
            print(f"[WARN] {year}: 未找到ALT文件")
            continue

        output_file = os.path.join(ALT_OUTPUT_DIR, f"ALT_{year}.tif")

        try:
            with rasterio.open(input_file) as src:
                # 确保CRS一致
                basin_reproj = basin.to_crs(src.crs)

                # 获取原始 nodata 值
                src_nodata = src.nodata if src.nodata is not None else -9999.0

                # 裁剪 - 保持原始 nodata 值
                out_image, out_transform = mask(
                    src,
                    basin_reproj.geometry,
                    crop=True,
                    nodata=src_nodata  # 使用原始的 nodata 值
                )

                # 更新元数据
                out_meta = src.meta.copy()
                out_meta.update({
                    "driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform,
                    "nodata": src_nodata,
                    "dtype": "float32"
                })

                # 保存
                with rasterio.open(output_file, "w", **out_meta) as dest:
                    dest.write(out_image)

            print(f"[OK] {year}: {output_file}")
            success_count += 1

        except Exception as e:
            print(f"[ERROR] {year}: 裁剪失败 - {e}")

    print(f"\n完成！成功裁剪 {success_count}/{len(YEARS)} 年")
    return success_count > 0


# ============================================================
# 可视化检查
# ============================================================
def visualize_alt(year=2010):
    """可视化某一年的ALT数据"""

    alt_file = os.path.join(ALT_OUTPUT_DIR, f"ALT_{year}.tif")

    if not os.path.exists(alt_file):
        print(f"文件不存在: {alt_file}")
        return

    with rasterio.open(alt_file) as src:
        alt_data = src.read(1)

    # 创建图形
    fig, ax = plt.subplots(figsize=(10, 8))

    # 绘制ALT
    im = ax.imshow(alt_data, cmap='viridis')
    plt.colorbar(im, ax=ax, label='Active Layer Thickness (m)')

    ax.set_title(f'目标流域活动层厚度 - {year}年')
    ax.set_xlabel('列')
    ax.set_ylabel('行')

    # 保存图片
    output_fig = os.path.join(PROJECT_ROOT, "results", "figures", f"ALT_{year}.png")
    os.makedirs(os.path.dirname(output_fig), exist_ok=True)
    plt.savefig(output_fig, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"图片已保存: {output_fig}")


# ============================================================
# 主函数
# ============================================================
if __name__ == "__main__":
    # 检查输入目录
    print("检查输入数据...")

    if not os.path.exists(BASIN_SHP):
        print(f"\n[WARN] 请先准备流域边界: {BASIN_SHP}")

    alt_files = [f for f in os.listdir(ALT_INPUT_DIR) if f.endswith('.tif')] if os.path.exists(ALT_INPUT_DIR) else []
    if not alt_files:
        print(f"\n[WARN] 请先将青藏高原ALT数据放入: {ALT_INPUT_DIR}")
    else:
        print(f"   找到 {len(alt_files)} 个ALT文件")

    # 如果数据都准备好了，执行裁剪
    if os.path.exists(BASIN_SHP) and alt_files:
        if clip_alt_to_basin():
            visualize_alt(2010)
    else:
        print("\n请准备好数据后再运行此脚本")

