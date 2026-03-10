# Reproducibility

This repository is designed for reproducible public workflows within a clearly defined scope.

## What The Repository Guarantees

- a stable public repository structure
- a documented config schema and template
- documented workflow entrypoints
- repository-level verification via `tools/verify_public_repo.py`
- a lightweight smoke test via `tools/smoke_test.py --lightweight`

## What The Repository Does Not Guarantee

- reproduction of private basin results without the corresponding private data
- publication of raw forcing archives or private run snapshots
- scientific equivalence across all new basins without basin-specific validation

## Reproducibility Workflow

Before publishing or sharing changes:

```bash
python tools/verify_public_repo.py
python tools/smoke_test.py --lightweight
```

These checks validate:

- public metadata and config shape
- absence of blocked repository terms
- availability of main CLI entrypoints
- integrity of the synthetic smoke-test example

## Recommended Handoff Artifact

For improved HBV runs, `metadata.json` is the main public machine-readable handoff artifact for downstream ALT-coupled workflows.

## Reproducibility Boundary

This is a method repository, not a private archive mirror. Reproducibility here means reproducible public interfaces and public workflow logic inside the repository boundary.
