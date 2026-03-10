# Permafrost ALT-Coupled HBV

Reusable Original / Improved / ALT-coupled HBV workflows for cold-region basin studies.

中文说明：这是一个面向公开复用的多年冻土 ALT 耦合 HBV 方法仓库，不包含特定流域案例叙述、私有数据或论文快照。

## What This Repository Is

This repository packages three related hydrological workflows:

1. `original_hbv`: a baseline distributed HBV configuration.
2. `improved_hbv`: a stronger baseline with zoned melt factors, glacier handling, source separation, and multi-objective calibration.
3. `alt_coupled_hbv`: a dynamic parameter-correction workflow driven by active layer thickness (ALT) reference fields and annual scale factors.

The public release is designed to be:

- reusable across basins
- configuration-driven
- free of basin-specific narrative material
- explicit about reproducibility boundaries
- easy to validate before publication

## What This Repository Is Not

- It is not a thesis archive.
- It does not ship raw forcing data, observed discharge, or run snapshots.
- It does not include basin-specific case-study text, figure captions, or default parameter sets tied to a private study area.
- It does not promise one-click execution for every external data source; some download steps require optional tools and external credentials.

## Model Tiers

### 1. Original HBV

Use this workflow for a structural baseline. It is intentionally simpler and serves as the main comparison target for the improved model.

### 2. Improved HBV

Use this workflow as the main calibration baseline before ALT coupling. It supports:

- zoned `CFMAX`
- glacier-aware runoff handling
- runoff source separation
- multi-objective calibration outputs such as `metadata.json`, `diagnostics.json`, and `parameters.txt`

### 3. ALT-Coupled HBV

Use this workflow after the improved HBV baseline is calibrated. It applies ALT-driven corrections to selected parameters and compares multiple coupling schemes using a shared output contract.

## Repository Layout

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
|- CONTRIBUTING.md
|- THIRD_PARTY.md
|- environment.yml
|- requirements-core.txt
`- LICENSE
```

## Installation

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

Optional tooling such as `whitebox`, `cdsapi`, `gdown`, or ArcGIS/ArcPy is documented in [docs/workflow_overview.md](docs/workflow_overview.md).

## Quick Start

1. Copy the basin template in `config/`.
2. Fill in basin paths, time ranges, observation settings, and ALT reference settings.
3. Run the data-preparation steps under `workflow/data_prep/`.
4. Calibrate the improved HBV baseline.
5. Run the original HBV baseline for comparison.
6. Generate ALT correction fields and run the ALT-coupled workflow.
7. Inspect standardized outputs and generic analysis scripts under `workflow/analysis/`.

Typical entrypoints:

```bash
python workflow/models/calibrate_improved_hbv.py --config config/basin.template.json
python workflow/models/run_original_hbv.py --config config/basin.template.json
python workflow/models/run_alt_coupled_hbv.py --config config/basin.template.json --improved-run path/to/metadata.json
```

## Validation Before You Push

Run these checks before publishing a change:

```bash
python tools/verify_public_repo.py
python tools/smoke_test.py --lightweight
```

For a short release checklist, see [docs/release_checklist.md](docs/release_checklist.md).

## Configuration Philosophy

The public release is explicit rather than heuristic:

- `observation_mode` is declared in the config
- seasonal objective windows are declared by month list
- ALT reference years are declared in the config
- coupling onset timing is declared in the config

This avoids hidden behavior tied to a single study basin.

## Standard Outputs

The repository keeps a stable public output contract:

- Original HBV: `original_hbv_summary.json`, `original_hbv_timeseries.csv`
- Improved HBV: `metadata.json`, `diagnostics.json`, `parameters.txt`
- ALT-coupled HBV: `summary.json`, `<scheme>_timeseries.csv`, `<scheme>_param_stats.csv`

## Synthetic Example

`examples/synthetic_basin/` is reserved for a minimal smoke-test dataset. It exists to validate interfaces and file contracts, not to demonstrate scientific performance.

## Documentation Index

- [docs/model_overview.md](docs/model_overview.md)
- [docs/data_contract.md](docs/data_contract.md)
- [docs/workflow_overview.md](docs/workflow_overview.md)
- [docs/output_guide.md](docs/output_guide.md)
- [docs/limitations.md](docs/limitations.md)
- [docs/privacy_scope.md](docs/privacy_scope.md)
- [docs/release_checklist.md](docs/release_checklist.md)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## Citation

See [CITATION.cff](CITATION.cff).

## License

This repository is released under the GNU General Public License v3.0. See [LICENSE](LICENSE).
