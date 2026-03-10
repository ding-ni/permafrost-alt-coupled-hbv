# FAQ

## Is this a full case-study archive?

No. This repository is a public method release. It does not contain private basin narratives, private run snapshots, or raw scientific datasets.

## Can I use this repository on a new basin?

Yes. That is the intended use case. Start with:

1. `config/basin.template.json`
2. `docs/workflow_overview.md`
3. `workflow/data_prep/`

## Why is there a synthetic example?

The synthetic example is a repository fixture for smoke testing. It exists to validate paths, interfaces, and output contracts. It is not a scientific benchmark.

## Why is `legacy_core/` still present?

The repository preserves a minimal vendored workflow core so that the public-facing wrappers remain stable and do not depend on a private workspace layout.

## Is the repository production-ready?

No. It is organized as research software with explicit publication boundaries and verification tooling, not as a production hydrological service.

## Why are some dependencies optional?

Some preprocessing paths depend on GIS or data-acquisition tools that are environment-specific, such as `whitebox`, `cdsapi`, `gdown`, or `arcpy`.

## What should I do before pushing changes?

Run:

```bash
python tools/verify_public_repo.py
python tools/smoke_test.py --lightweight
```

Then review [release_checklist.md](release_checklist.md).
