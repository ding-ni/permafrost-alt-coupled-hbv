# Third-Party Code and Dependency Notes

This repository is a public-facing workflow assembled from repository-native wrappers, declared Python dependencies, and a vendored minimal legacy core.

## Vendored Internal Legacy Core

The `legacy_core/` directory is not an external package mirror. It is a curated snapshot of workflow code carried into this repository to preserve the public interface and reproducibility boundary of this method release.

Current inventory:

1. `legacy_core/preprocess/*`
   Source: internal preprocessing workflow code from the original project workspace
   License in this repository: distributed under the repository GPL-3.0 release
   Reason retained: stable public preprocessing bridge for the wrapper entrypoints

2. `legacy_core/models/*`
   Source: internal HBV workflow code supporting original, improved, and ALT-coupled model execution
   License in this repository: distributed under the repository GPL-3.0 release
   Reason retained: preserve the public workflow contract without requiring a private upstream workspace layout

3. `legacy_core/glacier/*`
   Source: internal glacier-support scripts used by the public preprocessing workflow
   License in this repository: distributed under the repository GPL-3.0 release
   Reason retained: keep glacier preprocessing available to the public wrappers

## Declared Python Dependencies

Core Python dependencies are listed in:

- `environment.yml`
- `requirements-core.txt`

These are preferred over vendoring whenever package installation is sufficient.

## Optional External Tools

Some workflow steps may depend on optional tooling:

- `whitebox`
- `cdsapi`
- `gdown`
- ArcGIS / `arcpy`

These tools are optional because they are tied to specific GIS or data-acquisition environments.

## GPL Context

This public repository is released under GPL-3.0.

When adding external code or new vendored material:

- verify GPL compatibility
- preserve upstream notices where required
- document the source and reason for inclusion here

## Maintainer Rule

Do not vendor additional code without recording:

1. source project or origin
2. source URL if public
3. applicable license
4. snapshot or release identifier if available
5. why a declared dependency is not sufficient
