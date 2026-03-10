# Synthetic Basin Example

This directory contains a lightweight synthetic basin configuration for repository smoke checks.

It is intentionally minimal:

- it is not a scientific benchmark
- it does not represent a real basin
- it exists only to validate public interfaces, paths, and file contracts

Generate or refresh the example files with:

```bash
python tools/generate_synthetic_example.py --config examples/synthetic_basin/smoke_config.json
```
