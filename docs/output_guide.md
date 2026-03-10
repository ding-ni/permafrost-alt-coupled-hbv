# Output Guide

## Original HBV Outputs

These files support baseline comparison:

- `original_hbv_summary.json`
- `original_hbv_timeseries.csv`

## Improved HBV Outputs

These files support reproducibility and downstream ALT coupling:

- `metadata.json`
- `diagnostics.json`
- `parameters.txt`

`metadata.json` is the main machine-readable handoff artifact for later steps.

## ALT-Coupled Outputs

The ALT workflow produces:

- a workflow-level `summary.json`
- one timeseries file per scheme
- one parameter-statistics file per scheme when available

## Interpretation Rule

The public repository distinguishes:

- the objective-window metrics actually used in calibration
- the full-period metrics used for cross-scheme comparison

This distinction should remain visible in outputs and documentation.
