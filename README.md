# Permafrost ALT-Coupled HBV

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB.svg)](environment.yml)
[![CI](https://github.com/ding-ni/permafrost-alt-coupled-hbv/actions/workflows/ci.yml/badge.svg)](https://github.com/ding-ni/permafrost-alt-coupled-hbv/actions/workflows/ci.yml)
[![Status](https://img.shields.io/badge/status-research%20software-6f42c1.svg)](docs/README.md)

[English](README.md) | [简体中文](README.zh-CN.md)

Reusable Original / Improved / ALT-coupled HBV workflows for cold-region basin studies.

## Project Status

This repository is a public method release and research-software codebase.

- It is intended to be reusable across basins.
- It is actively curated as a public-facing workflow repository.
- It is not positioned as a production hydrological platform or a thesis archive.

Current repository signals:

- project status: [PROJECT_STATUS.md](PROJECT_STATUS.md)
- planned improvements: [ROADMAP.md](ROADMAP.md)
- release history: [CHANGELOG.md](CHANGELOG.md)

## Who This Repository Is For

This repository is primarily for:

- researchers working on cold-region hydrology
- students or collaborators reproducing the workflow on a new basin
- maintainers who need a basin-agnostic public version of the method

## Workflow At A Glance

The repository exposes three related model tiers:

1. `original_hbv`
This is the structural baseline used for comparison.

2. `improved_hbv`
This is the main calibration baseline with zoned melt factors, glacier-aware handling, source separation, and reproducible calibration outputs.

3. `alt_coupled_hbv`
This applies ALT-driven parameter corrections on top of the improved HBV baseline and compares multiple coupling schemes under a shared output contract.

## Start Here

If this is your first time in the repository, read these files in order:

1. [docs/README.md](docs/README.md)
2. [docs/workflow_overview.md](docs/workflow_overview.md)
3. [config/basin.template.json](config/basin.template.json)

If you are evaluating whether this repository fits your use case, read:

1. [PROJECT_STATUS.md](PROJECT_STATUS.md)
2. [docs/privacy_scope.md](docs/privacy_scope.md)
3. [docs/limitations.md](docs/limitations.md)

## Quick Start

### Installation

Conda:

```bash
conda env create -f environment.yml
conda activate permafrost-alt-hbv
```

pip:

```bash
python -m venv .venv
. .venv/Scripts/activate
pip install -r requirements-core.txt
```

### Minimal Working Sequence

1. Copy the basin config template in `config/`.
2. Fill in basin paths, time windows, observation settings, and ALT reference settings.
3. Run the data-preparation steps under `workflow/data_prep/`.
4. Calibrate the improved HBV baseline.
5. Run the original HBV baseline for comparison.
6. Generate ALT correction fields and run the ALT-coupled workflow.
7. Inspect outputs under `workflow/analysis/`.

Typical entrypoints:

```bash
python workflow/models/calibrate_improved_hbv.py --config config/basin.template.json
python workflow/models/run_original_hbv.py --config config/basin.template.json
python workflow/models/run_alt_coupled_hbv.py --config config/basin.template.json --improved-run path/to/metadata.json
```

## Public Scope And Reproducibility Boundary

This repository is designed to be:

- configuration-driven
- basin-agnostic
- explicit about public vs private material

This repository does not ship:

- raw forcing archives
- private observed discharge files
- private run snapshots
- thesis-only case-study prose
- case-specific figure chains tied to a private basin

For the exact publication boundary, see [docs/privacy_scope.md](docs/privacy_scope.md).

## Repository Map

```text
permafrost-alt-coupled-hbv/
|- config/       public configuration template and schema
|- docs/         method, workflow, scope, and release docs
|- examples/     synthetic smoke-test example only
|- legacy_core/  minimal vendored workflow core used by public wrappers
|- shared/       shared runtime helpers and alias logic
|- tools/        verification and smoke-test utilities
|- workflow/     public-facing preprocessing, model, and analysis entrypoints
`- .github/      issue and pull request templates
```

## Standard Output Contract

The repository keeps a stable public output contract:

- Original HBV: `original_hbv_summary.json`, `original_hbv_timeseries.csv`
- Improved HBV: `metadata.json`, `diagnostics.json`, `parameters.txt`
- ALT-coupled HBV: `summary.json`, `<scheme>_timeseries.csv`, `<scheme>_param_stats.csv`

## Validation Before You Push

Run these checks before publishing an update:

```bash
python tools/verify_public_repo.py
python tools/smoke_test.py --lightweight
```

For a release checklist, see [docs/release_checklist.md](docs/release_checklist.md).

## Documentation Index

- [docs/README.md](docs/README.md)
- [docs/faq.md](docs/faq.md)
- [docs/architecture.md](docs/architecture.md)
- [docs/reproducibility.md](docs/reproducibility.md)
- [docs/model_overview.md](docs/model_overview.md)
- [docs/data_contract.md](docs/data_contract.md)
- [docs/workflow_overview.md](docs/workflow_overview.md)
- [docs/output_guide.md](docs/output_guide.md)
- [docs/limitations.md](docs/limitations.md)
- [docs/privacy_scope.md](docs/privacy_scope.md)
- [docs/release_checklist.md](docs/release_checklist.md)

## Support

- Workflow questions: see [SUPPORT.md](SUPPORT.md)
- Contribution process: see [CONTRIBUTING.md](CONTRIBUTING.md)
- Security reporting guidance: see [SECURITY.md](SECURITY.md)
- Community expectations: see [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)

## Citation

See [CITATION.cff](CITATION.cff).

## License

This repository is released under the GNU General Public License v3.0. See [LICENSE](LICENSE).
