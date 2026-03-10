# Permafrost ALT-Coupled HBV

[English](README.md) | [简体中文](README.zh-CN.md)

一个面向寒区流域研究的 Original / Improved / ALT-coupled HBV 可复用方法仓库。

## 这个仓库是什么

这个仓库整理并公开了三套互相关联的水文工作流：

1. `original_hbv`：原始分布式 HBV 基线模型。
2. `improved_hbv`：增强版基线模型，加入分区融雪因子、冰川处理、径流来源分离和多目标率定。
3. `alt_coupled_hbv`：基于活动层厚度（ALT）参考场与年际尺度因子的动态参数修正工作流。

公开版的目标是：

- 可迁移到不同流域
- 由配置文件驱动
- 不包含特定流域的案例叙述
- 明确可复现边界
- 在公开发布前可以快速自检

## 这个仓库不是什么

- 不是论文归档目录
- 不提供原始气象强迫、实测径流或历史运行快照
- 不包含私有研究区的正文描述、图注或默认参数集
- 不保证所有外部数据源都能一键下载，因为其中一些步骤依赖可选工具和外部账号权限

## 三套模型层级

### 1. Original HBV

这是结构性基线模型，主要用于和改进型 HBV 做对照。

### 2. Improved HBV

这是 ALT 耦合前的主要率定基线，支持：

- 分区 `CFMAX`
- 冰川相关产流处理
- 径流来源分离
- 标准率定输出，例如 `metadata.json`、`diagnostics.json`、`parameters.txt`

### 3. ALT-Coupled HBV

在改进型 HBV 完成率定后运行。它会对选定参数进行 ALT 驱动修正，并在统一输出契约下比较多种耦合方案。

## 仓库结构

```text
permafrost-alt-coupled-hbv/
|- config/
|- docs/
|- examples/
|- legacy_core/
|- shared/
|- tools/
|- workflow/
|- .github/
|- README.md
|- README.zh-CN.md
|- CONTRIBUTING.md
|- THIRD_PARTY.md
|- environment.yml
|- requirements-core.txt
`- LICENSE
```

## 安装

### Conda

```bash
conda env create -f environment.yml
conda activate permafrost-alt-hbv
```

### pip

```bash
python -m venv .venv
. .venv/Scripts/activate
pip install -r requirements-core.txt
```

`whitebox`、`cdsapi`、`gdown`、ArcGIS/ArcPy 等可选依赖见 [docs/workflow_overview.md](docs/workflow_overview.md)。

## 快速开始

1. 复制 `config/` 下的流域模板配置。
2. 填写流域路径、时间范围、观测口径和 ALT 参考期。
3. 按顺序运行 `workflow/data_prep/` 下的数据准备脚本。
4. 率定改进型 HBV。
5. 运行原始 HBV 作为对照。
6. 生成 ALT 修正参数并运行 ALT 耦合工作流。
7. 使用 `workflow/analysis/` 下的通用分析脚本检查输出。

典型入口：

```bash
python workflow/models/calibrate_improved_hbv.py --config config/basin.template.json
python workflow/models/run_original_hbv.py --config config/basin.template.json
python workflow/models/run_alt_coupled_hbv.py --config config/basin.template.json --improved-run path/to/metadata.json
```

## 发布前自检

建议在每次公开更新前运行：

```bash
python tools/verify_public_repo.py
python tools/smoke_test.py --lightweight
```

简短发布清单见 [docs/release_checklist.md](docs/release_checklist.md)。

## 配置理念

公开版尽量避免隐式启发式规则，而是显式声明：

- `observation_mode`
- 目标窗口月份
- ALT 参考年份
- 耦合起始年份

这样可以避免把单一流域的隐含假设带进公开方法仓库。

## 标准输出

仓库当前维护以下稳定输出契约：

- Original HBV：`original_hbv_summary.json`、`original_hbv_timeseries.csv`
- Improved HBV：`metadata.json`、`diagnostics.json`、`parameters.txt`
- ALT-coupled HBV：`summary.json`、`<scheme>_timeseries.csv`、`<scheme>_param_stats.csv`

## 合成示例

`examples/synthetic_basin/` 是一个最小化 smoke-test 示例，只用于验证接口和文件契约，不代表真实科学结果。

## 文档索引

- [docs/model_overview.md](docs/model_overview.md)
- [docs/data_contract.md](docs/data_contract.md)
- [docs/workflow_overview.md](docs/workflow_overview.md)
- [docs/output_guide.md](docs/output_guide.md)
- [docs/limitations.md](docs/limitations.md)
- [docs/privacy_scope.md](docs/privacy_scope.md)
- [docs/release_checklist.md](docs/release_checklist.md)

## 贡献

见 [CONTRIBUTING.md](CONTRIBUTING.md)。

## 引用

见 [CITATION.cff](CITATION.cff)。

## 许可证

本仓库采用 GNU General Public License v3.0。见 [LICENSE](LICENSE)。
