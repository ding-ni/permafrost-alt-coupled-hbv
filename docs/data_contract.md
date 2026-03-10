# Data Contract

This document defines the public input and output expectations for the reusable workflow.

## Required Inputs

### Basin geometry

- basin boundary shapefile or equivalent vector dataset
- DEM raster covering the full basin

### Hydrometeorological data

- observed discharge CSV
- temperature forcing
- precipitation forcing
- potential evapotranspiration inputs or prepared PET rasters

### Cryosphere data

- ALT rasters for the target years
- optional glacier boundary vector data

## Configuration Contract

The authoritative public config artifacts are:

- `config/basin.template.json`
- `config/basin.schema.json`

At minimum, the public config is expected to declare:

- basin identity
- workspace root
- time windows
- observation mode
- objective months when seasonal evaluation is used
- forcing preferences
- ALT reference period
- coupling start year
- improved HBV handoff path when ALT-coupled workflows are used

## Output Contract

### Original HBV

- `original_hbv_summary.json`
- `original_hbv_timeseries.csv`

### Improved HBV

- `metadata.json`
- `diagnostics.json`
- `parameters.txt`

### ALT-coupled HBV

- `summary.json`
- `<scheme>_timeseries.csv`
- `<scheme>_param_stats.csv`

## Public Repository Rule

The repository stores interface definitions and fixtures, not real scientific datasets. Do not commit:

- raw forcing archives
- private observation files
- run snapshots
- site-specific figures
