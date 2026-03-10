# -*- coding: utf-8 -*-
"""
阶段1：GIS 水文分析
生成流向栅格和流量累积栅格

支持两种方式：
1. WhiteboxTools (推荐)
2. ArcPy (备选)

输入：
- data/gis/dem_1km.tif  DEM栅格

输出：
- data/gis/dem_filled.tif      填洼后的DEM
- data/gis/flow_direction.tif  D8流向栅格
- data/gis/flow_accumulation.tif 流量累积栅格
"""

import sys
import os

# 添加 Hapi 路径
sys.path.insert(0, r".\Hapi-main\Hapi-main\src")

# ============================================================
# 路径配置
# ============================================================
GIS_DIR = r".\data\gis"

# 输入文件
DEM_INPUT = os.path.join(GIS_DIR, "dem_1km.tif")

# 输出文件
DEM_FILLED = os.path.join(GIS_DIR, "dem_filled.tif")
FLOW_DIR = os.path.join(GIS_DIR, "flow_direction.tif")
FLOW_ACC = os.path.join(GIS_DIR, "flow_accumulation.tif")


# ============================================================
# 方法1：使用 WhiteboxTools
# ============================================================
def run_whitebox():
    """使用 WhiteboxTools 进行水文分析"""
    from whitebox import WhiteboxTools

    wbt = WhiteboxTools()
    wbt.work_dir = GIS_DIR
    wbt.verbose = True

    print("=" * 60)
    print("使用 WhiteboxTools 进行水文分析")
    print("=" * 60)

    # 检查输入文件
    if not os.path.exists(DEM_INPUT):
        print(f"[ERROR] DEM文件不存在: {DEM_INPUT}")
        print("   请先将DEM放入 data/gis/ 目录")
        return False

    # 步骤1：填充洼地
    print("\n[1/3] 填充洼地...")
    wbt.fill_depressions(
        dem=os.path.basename(DEM_INPUT),
        output=os.path.basename(DEM_FILLED)
    )
    print(f"   输出: {DEM_FILLED}")

    # 步骤2：生成 D8 流向
    print("\n[2/3] 生成 D8 流向...")
    wbt.d8_pointer(
        dem=os.path.basename(DEM_FILLED),
        output=os.path.basename(FLOW_DIR)
    )
    print(f"   输出: {FLOW_DIR}")

    # 步骤3：生成流量累积
    print("\n[3/3] 生成流量累积...")
    wbt.d8_flow_accumulation(
        i=os.path.basename(DEM_FILLED),
        output=os.path.basename(FLOW_ACC),
        out_type="cells"  # 输出单元格数量
    )
    print(f"   输出: {FLOW_ACC}")

    print("\n[OK] WhiteboxTools 水文分析完成！")
    return True


# ============================================================
# 方法2：使用 ArcPy
# ============================================================
def run_arcpy():
    """使用 ArcPy 进行水文分析"""
    import arcpy
    from arcpy.sa import Fill, FlowDirection, FlowAccumulation

    # 启用空间分析扩展
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = GIS_DIR

    print("=" * 60)
    print("使用 ArcPy 进行水文分析")
    print("=" * 60)

    # 检查输入文件
    if not os.path.exists(DEM_INPUT):
        print(f"[ERROR] DEM文件不存在: {DEM_INPUT}")
        print("   请先将DEM放入 data/gis/ 目录")
        return False

    # 步骤1：填充洼地
    print("\n[1/3] 填充洼地 (Fill)...")
    dem_filled = Fill(DEM_INPUT)
    dem_filled.save(DEM_FILLED)
    print(f"   输出: {DEM_FILLED}")

    # 步骤2：生成流向
    print("\n[2/3] 生成流向 (D8)...")
    flow_dir = FlowDirection(DEM_FILLED, "NORMAL")
    flow_dir.save(FLOW_DIR)
    print(f"   输出: {FLOW_DIR}")

    # 步骤3：生成流量累积
    print("\n[3/3] 生成流量累积...")
    flow_acc = FlowAccumulation(FLOW_DIR)
    flow_acc.save(FLOW_ACC)
    print(f"   输出: {FLOW_ACC}")

    # 归还扩展许可
    arcpy.CheckInExtension("Spatial")

    print("\n[OK] ArcPy 水文分析完成！")
    return True


# ============================================================
# 验证结果
# ============================================================
def verify_results():
    """验证生成的栅格"""
    import rasterio

    print("\n" + "=" * 60)
    print("验证结果")
    print("=" * 60)

    files = [
        ("DEM (填洼后)", DEM_FILLED),
        ("流向栅格", FLOW_DIR),
        ("流量累积", FLOW_ACC),
    ]

    for name, path in files:
        if os.path.exists(path):
            with rasterio.open(path) as src:
                print(f"\n[OK] {name}")
                print(f"   路径: {path}")
                print(f"   尺寸: {src.width} x {src.height}")
                print(f"   CRS: {src.crs}")
                print(f"   分辨率: {src.res}")
        else:
            print(f"\n[ERROR] {name} - 文件不存在")


# ============================================================
# 主函数
# ============================================================
if __name__ == "__main__":
    print("目标流域水文分析")
    print("=" * 60)

    # 检查DEM是否存在
    if not os.path.exists(DEM_INPUT):
        print(f"\n[WARN] DEM文件尚未准备: {DEM_INPUT}")
        print("请将目标流域DEM放入 data/gis/ 目录，命名为 dem_1km.tif")
        print("\n脚本将在DEM准备好后运行")
        sys.exit(0)

    # 尝试使用 WhiteboxTools
    try:
        from whitebox import WhiteboxTools
        success = run_whitebox()
    except ImportError:
        print("WhiteboxTools 未安装，尝试使用 ArcPy...")
        try:
            import arcpy
            success = run_arcpy()
        except ImportError:
            print("[ERROR] 需要安装 WhiteboxTools 或 ArcPy")
            print("   pip install whitebox")
            sys.exit(1)

    # 验证结果
    if success:
        verify_results()

