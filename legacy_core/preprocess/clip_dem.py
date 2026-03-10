# -*- coding: utf-8 -*-
"""
阶段0B：裁剪青藏高原 DEM 到目标流域

从青藏高原整体DEM裁剪出流域范围的DEM

输入：
- E:/permaforst/zuankong/dem_best/elevation.tif  青藏高原DEM
- data/gis/basin.shp                     流域边界

输出：
- data/gis/dem_1km.tif                            裁剪后的流域DEM
"""

import sys
import os

sys.path.insert(0, r".\Hapi-main\Hapi-main\src")

import numpy as np
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
import geopandas as gpd

# ============================================================
# 路径配置
# ============================================================
PROJECT_ROOT = r"."

# 输入
TP_DEM_FILE = r"E:\permaforst\zuankong\dem_best\elevation.tif"
BASIN_SHP = os.path.join(PROJECT_ROOT, "data", "gis", "basin.shp")

# 输出
OUTPUT_DEM = os.path.join(PROJECT_ROOT, "data", "gis", "dem_1km.tif")


# ============================================================
# 裁剪 DEM
# ============================================================
def clip_dem_to_basin():
    """裁剪青藏高原DEM到流域范围"""

    print("=" * 60)
    print("裁剪青藏高原 DEM 到目标流域")
    print("=" * 60)

    # 检查输入文件
    if not os.path.exists(TP_DEM_FILE):
        print(f"\n[ERROR] 青藏高原DEM不存在: {TP_DEM_FILE}")
        return False

    if not os.path.exists(BASIN_SHP):
        print(f"\n[ERROR] 流域边界不存在: {BASIN_SHP}")
        print("请先准备目标流域边界 Shapefile")
        return False

    # 读取流域边界
    print(f"\n[1/3] 读取流域边界...")
    basin = gpd.read_file(BASIN_SHP)
    print(f"   流域面积: {basin.geometry.area.sum() / 1e6:.2f} km2")
    print(f"   CRS: {basin.crs}")

    # 打开青藏高原DEM
    print(f"\n[2/3] 读取青藏高原DEM...")
    with rasterio.open(TP_DEM_FILE) as src:
        print(f"   原始DEM尺寸: {src.width} x {src.height}")
        print(f"   原始分辨率: {src.res}")
        print(f"   原始CRS: {src.crs}")

        # 确保流域边界与DEM CRS一致
        basin_reproj = basin.to_crs(src.crs)

        # 裁剪
        print(f"\n[3/3] 裁剪DEM到流域范围...")
        out_image, out_transform = mask(
            src,
            basin_reproj.geometry,
            crop=True,
            nodata=-9999  # 使用 -9999 作为 nodata（适合整数类型）
        )

        # 更新元数据
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
            "nodata": -9999,
            "compress": "lzw"  # 压缩以节省空间
        })

        # 保存
        os.makedirs(os.path.dirname(OUTPUT_DEM), exist_ok=True)

        with rasterio.open(OUTPUT_DEM, "w", **out_meta) as dst:
            dst.write(out_image)

    print(f"\n[OK] DEM已保存: {OUTPUT_DEM}")

    # 验证结果
    with rasterio.open(OUTPUT_DEM) as dem:
        data = dem.read(1)
        print(f"\n裁剪后DEM信息:")
        print(f"   尺寸: {dem.width} x {dem.height}")
        print(f"   分辨率: {dem.res}")
        print(f"   CRS: {dem.crs}")
        print(f"   高程范围: {np.nanmin(data):.1f} ~ {np.nanmax(data):.1f} m")
        print(f"   平均高程: {np.nanmean(data):.1f} m")

    return True


# ============================================================
# 重采样到 1km (如果需要)
# ============================================================
def resample_to_1km():
    """
    将DEM重采样到 1km 分辨率
    注意：只有当原始分辨率不是1km时才需要运行
    """

    print("\n" + "=" * 60)
    print("重采样 DEM 到 1km 分辨率")
    print("=" * 60)

    if not os.path.exists(OUTPUT_DEM):
        print("[ERROR] 请先运行 clip_dem_to_basin()")
        return False

    # 读取原始DEM
    with rasterio.open(OUTPUT_DEM) as src:
        current_res = src.res[0]
        print(f"当前分辨率: {current_res} 度")

        # 检查是否需要重采样
        # 1km ≈ 0.009度（在青藏高原纬度约34°）
        target_res = 0.009  # 度
        tolerance = 0.001

        if abs(current_res - target_res) < tolerance:
            print(f"分辨率已经接近1km ({current_res:.6f}度)，无需重采样")
            return True

        print(f"需要重采样: {current_res:.6f}度 -> {target_res:.6f}度")

        # 计算新的transform
        transform, width, height = calculate_default_transform(
            src.crs, src.crs,
            src.width, src.height,
            *src.bounds,
            resolution=(target_res, target_res)
        )

        # 更新元数据
        kwargs = src.meta.copy()
        kwargs.update({
            'transform': transform,
            'width': width,
            'height': height
        })

        # 重采样
        print("重采样中...")
        temp_file = OUTPUT_DEM.replace(".tif", "_1km.tif")

        with rasterio.open(temp_file, 'w', **kwargs) as dst:
            for band in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, band),
                    destination=rasterio.band(dst, band),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=src.crs,
                    resampling=Resampling.bilinear
                )

        # 替换原文件
        os.replace(temp_file, OUTPUT_DEM)
        print(f"[OK] 重采样完成: {width} x {height}")

    return True


# ============================================================
# 可视化检查
# ============================================================
def visualize_dem():
    """可视化DEM"""

    if not os.path.exists(OUTPUT_DEM):
        print("DEM文件不存在，跳过可视化")
        return

    import matplotlib.pyplot as plt

    with rasterio.open(OUTPUT_DEM) as src:
        dem_data = src.read(1)

    # 创建图形
    fig, ax = plt.subplots(figsize=(10, 8))

    im = ax.imshow(dem_data, cmap='terrain')
    plt.colorbar(im, ax=ax, label='Elevation (m)')

    ax.set_title('Example Basin DEM')
    ax.set_xlabel('Column')
    ax.set_ylabel('Row')

    # 保存
    output_fig = os.path.join(PROJECT_ROOT, "results", "figures", "dem_basin.png")
    os.makedirs(os.path.dirname(output_fig), exist_ok=True)
    plt.savefig(output_fig, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\n[OK] DEM可视化已保存: {output_fig}")


# ============================================================
# 主函数
# ============================================================
if __name__ == "__main__":
    print("目标流域 DEM 裁剪工具")
    print("=" * 60)

    # 检查源文件
    print("\n检查源文件:")
    print(f"   青藏高原DEM: {TP_DEM_FILE}")
    print(f"   存在: {os.path.exists(TP_DEM_FILE)}")
    print(f"\n   流域边界: {BASIN_SHP}")
    print(f"   存在: {os.path.exists(BASIN_SHP)}")

    if not os.path.exists(TP_DEM_FILE) or not os.path.exists(BASIN_SHP):
        print("\n[WARNING] 缺少必要文件，请准备后再运行")
        sys.exit(0)

    # 执行裁剪
    success = clip_dem_to_basin()

    if success:
        # 重采样（如需要）
        resample_to_1km()

        # 可视化
        visualize_dem()

        print("\n" + "=" * 60)
        print("[DONE] DEM处理完成！")
        print("=" * 60)
        print(f"\n输出文件: {OUTPUT_DEM}")

