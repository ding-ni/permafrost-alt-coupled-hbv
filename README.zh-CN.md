# Permafrost ALT-Coupled HBV

[English](README.md) | [简体中文](README.zh-CN.md)

一个面向寒区流域研究的 Original / Improved / ALT-coupled HBV 可复用方法仓库。

## 项目状态

这个仓库当前定位为“公开方法发布 + 研究型软件仓库”。

- 面向新流域复用，而不是单一案例存档
- 正在按公开工程仓库的标准持续整理
- 不把自己包装成生产级水文平台，也不是论文归档目录

当前仓库关键信号：

- 项目状态： [PROJECT_STATUS.md](PROJECT_STATUS.md)
- 后续计划： [ROADMAP.md](ROADMAP.md)
- 更新记录： [CHANGELOG.md](CHANGELOG.md)

## 适用对象

这个仓库主要适合：

- 寒区水文研究人员
- 需要把方法迁移到新流域的学生或合作作者
- 需要维护一个不依赖具体案例的公开方法仓库的维护者

## 工作流总览

仓库公开了三层相关模型：

1. `original_hbv`
结构性基线模型，用于对照。

2. `improved_hbv`
主要率定基线，包含分区融雪因子、冰川处理、径流来源分离和可复现率定输出。

3. `alt_coupled_hbv`
在改进型 HBV 基础上做 ALT 驱动参数修正，并比较多种耦合方案。

## 从哪里开始读

如果你第一次进入这个仓库，建议按这个顺序看：

1. [docs/README.md](docs/README.md)
2. [docs/workflow_overview.md](docs/workflow_overview.md)
3. [config/basin.template.json](config/basin.template.json)

如果你想先判断这个仓库是否适合你的使用场景，建议先看：

1. [PROJECT_STATUS.md](PROJECT_STATUS.md)
2. [docs/privacy_scope.md](docs/privacy_scope.md)
3. [docs/limitations.md](docs/limitations.md)

## 快速开始

### 安装

Conda：

```bash
conda env create -f environment.yml
conda activate permafrost-alt-hbv
```

pip：

```bash
python -m venv .venv
. .venv/Scripts/activate
pip install -r requirements-core.txt
```

### 最小工作顺序

1. 复制 `config/` 下的流域配置模板。
2. 填写流域路径、时间窗口、观测口径和 ALT 参考期。
3. 运行 `workflow/data_prep/` 下的数据准备脚本。
4. 率定改进型 HBV。
5. 运行原始 HBV 作为对照。
6. 生成 ALT 修正参数并运行 ALT 耦合工作流。
7. 使用 `workflow/analysis/` 检查输出。

典型入口：

```bash
python workflow/models/calibrate_improved_hbv.py --config config/basin.template.json
python workflow/models/run_original_hbv.py --config config/basin.template.json
python workflow/models/run_alt_coupled_hbv.py --config config/basin.template.json --improved-run path/to/metadata.json
```

## 公开范围与可复现边界

这个仓库的原则是：

- 配置驱动
- 不绑定具体流域
- 明确区分公开内容和私有内容

这个仓库不提供：

- 原始气象强迫档案
- 私有实测径流文件
- 私有运行快照
- 论文专用案例文字
- 绑定私有流域的图件主链

详细边界见 [docs/privacy_scope.md](docs/privacy_scope.md)。

## 仓库结构

```text
permafrost-alt-coupled-hbv/
|- config/       公开配置模板与 schema
|- docs/         方法、流程、边界与发布文档
|- examples/     仅用于 smoke-test 的合成示例
|- legacy_core/  为公开包装层保留的最小核心代码
|- shared/       共享运行工具与公开别名逻辑
|- tools/        仓库校验与 smoke-test 工具
|- workflow/     公开入口脚本
`- .github/      Issue / PR 模板
```

## 标准输出契约

仓库当前维护以下稳定输出：

- Original HBV：`original_hbv_summary.json`、`original_hbv_timeseries.csv`
- Improved HBV：`metadata.json`、`diagnostics.json`、`parameters.txt`
- ALT-coupled HBV：`summary.json`、`<scheme>_timeseries.csv`、`<scheme>_param_stats.csv`

## 发布前校验

在公开更新前建议运行：

```bash
python tools/verify_public_repo.py
python tools/smoke_test.py --lightweight
```

发布检查清单见 [docs/release_checklist.md](docs/release_checklist.md)。

## 文档索引

- [docs/README.md](docs/README.md)
- [docs/faq.md](docs/faq.md)
- [docs/model_overview.md](docs/model_overview.md)
- [docs/data_contract.md](docs/data_contract.md)
- [docs/workflow_overview.md](docs/workflow_overview.md)
- [docs/output_guide.md](docs/output_guide.md)
- [docs/limitations.md](docs/limitations.md)
- [docs/privacy_scope.md](docs/privacy_scope.md)
- [docs/release_checklist.md](docs/release_checklist.md)

## 支持与维护

- 工作流使用问题：见 [SUPPORT.md](SUPPORT.md)
- 贡献方式：见 [CONTRIBUTING.md](CONTRIBUTING.md)
- 安全问题上报：见 [SECURITY.md](SECURITY.md)
- 社区行为规范：见 [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)

## 引用

见 [CITATION.cff](CITATION.cff)。

## 许可证

本仓库采用 GNU General Public License v3.0。见 [LICENSE](LICENSE)。
