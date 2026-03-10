# Workflow Overview

## End-to-End Sequence

1. Prepare basin geometry and workspace paths.
2. Prepare DEM, flow routing, and elevation products.
3. Prepare meteorological forcing rasters.
4. Prepare glacier products when applicable.
5. Prepare ALT reference fields and correction products.
6. Run the improved HBV calibration baseline.
7. Run the original HBV baseline for comparison.
8. Run ALT-coupled multi-scheme analysis.
9. Inspect outputs with the generic analysis tools.

## External Tooling

Some data-preparation steps may require optional tools:

- `whitebox` for routing preprocessing
- ArcGIS / `arcpy` as an alternative GIS backend
- `cdsapi` for CDS-backed downloads
- `gdown` for Google Drive based file retrieval

These tools are optional because not every basin workflow will use every acquisition path.

## Recommended Practice

- Use explicit config fields instead of implicit heuristics.
- Keep all basin-specific data outside the repository.
- Use a separate workspace per basin.
- Treat the synthetic example as an interface test, not a scientific benchmark.
- Run `python tools/verify_public_repo.py` and `python tools/smoke_test.py --lightweight` before pushing public updates.

## Release Boundary

The public workflow is method-oriented. If you need a paper-ready figure chain for a specific basin, maintain it outside this repository or rebuild it using the generic analysis outputs.
