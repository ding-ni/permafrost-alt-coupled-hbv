# Contributing

Thank you for improving this repository.

This repository is maintained as a public, basin-agnostic method release. Contributions should strengthen reuse, clarity, and reproducibility without reintroducing private case-study material.

## Ground Rules

- Keep the repository basin-agnostic.
- Do not add private raw data, private run snapshots, or private case-study text.
- Keep public interfaces and CLI surfaces in English.
- Chinese documentation support is allowed when it improves accessibility for the intended audience.
- Preserve the public output contracts unless there is a documented migration reason.

## Good Contribution Areas

- documentation improvements
- workflow bug fixes
- configuration/schema improvements
- basin-agnostic analysis utilities
- reproducibility and verification tooling

## Before You Open A Pull Request

1. Run `python tools/verify_public_repo.py`
2. Run `python tools/smoke_test.py --lightweight`
3. Search for private basin identifiers, local absolute paths, and private filenames
4. Update the relevant docs if workflow behavior, config schema, or outputs changed

## Documentation Changes

When changing docs:

- prefer concise, operational language
- state prerequisites and limitations clearly
- keep English primary for interfaces and top-level engineering docs
- keep Chinese support concise and intentional

## Code Changes

When changing workflow code:

- prefer explicit config over hidden heuristics
- avoid hardcoded absolute paths
- keep optional dependencies optional
- document new CLI arguments in `README.md` or `docs/`

## Pull Request Checklist

- [ ] No private data or private run outputs added
- [ ] No basin-specific content added
- [ ] Docs updated if public behavior changed
- [ ] Output contracts preserved or explicitly documented
- [ ] New dependencies documented
