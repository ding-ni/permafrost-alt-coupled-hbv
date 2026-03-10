# Model Overview

This repository exposes three related workflows that share a common data contract and workspace layout.

简述：三套模型是递进关系，不是互相独立的三个项目。

## Original HBV

Purpose:

- provide a simple structural baseline
- support model-comparison reporting

Typical use:

- run after data preparation
- compare against the improved HBV baseline

## Improved HBV

Purpose:

- serve as the main calibrated baseline for reusable applications
- provide the parameter set required by ALT coupling

Expected capabilities:

- zoned melt-factor behavior
- glacier-aware processing
- runoff component tracing
- reproducible calibration metadata

## ALT-Coupled HBV

Purpose:

- evaluate dynamic parameter corrections linked to active layer thickness
- compare multiple coupling schemes under a single workflow

Expected inputs:

- improved HBV calibration result
- ALT reference field and annual scaling products

## Shared Design Principles

- configuration-driven paths
- basin-agnostic documentation
- explicit reproducibility boundaries
- stable output filenames

## Public Release Boundary

This documentation explains the method family, not a single historical case study. The public repository intentionally excludes:

- basin-specific prose
- thesis-only figures
- private run directories
- private forcing and discharge datasets
