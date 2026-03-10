# Release Checklist

Use this checklist before publishing a major update.

## Repository Hygiene

- Run `python tools/verify_public_repo.py`
- Run `python tools/smoke_test.py --lightweight`
- Confirm `git status` is clean before tagging or pushing
- Confirm no raw basin data, observed discharge archives, run snapshots, or private notes were added

## Public Scope

- Confirm no basin-specific narrative text was introduced
- Confirm no local absolute paths were introduced
- Confirm no private credentials or API config files are tracked
- Confirm new output examples are synthetic or intentionally public

## Documentation

- Update `README.md` if the workflow, required inputs, or outputs changed
- Update files in `docs/` if public interfaces changed
- Update `CITATION.cff` if authorship or citation text changed
- Update `THIRD_PARTY.md` if vendored or external dependencies changed

## Reproducibility

- Confirm the basin config schema still matches the public config template
- Confirm main CLI entrypoints still expose `--help`
- Confirm the synthetic example still runs through the lightweight smoke checks

## Recommended Push Flow

```bash
git add .
git commit -m "Describe the release change"
git push
```
