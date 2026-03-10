# Third-Party Code and Dependency Notes

This repository packages a public-facing workflow built from maintained project code and a minimal legacy core.

## Current Policy

- Keep third-party code use explicit.
- Prefer declared package dependencies over vendoring.
- Vendor code only when it is necessary for reproducibility or interface stability.
- Record the source, license, and reason for inclusion whenever vendored code is added.

## Expected Dependency Categories

### Core runtime packages

Core Python dependencies are listed in:

- `environment.yml`
- `requirements-core.txt`

### Optional external tools

Some workflow steps may depend on optional tooling:

- `whitebox`
- `cdsapi`
- `gdown`
- ArcGIS / `arcpy`

These are optional because they are tied to specific data-acquisition or GIS environments.

## GPL Context

This public repository is released under GPL-3.0.

When adding third-party code:

- verify license compatibility with GPL-3.0
- preserve upstream notices where required
- document any local modifications

## Maintainer Checklist

Before adding any vendored code, record:

1. upstream project name
2. upstream source URL
3. upstream license
4. commit, release, or snapshot identifier
5. why package dependency alone is not sufficient
