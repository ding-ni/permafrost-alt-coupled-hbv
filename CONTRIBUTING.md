# Contributing

Thank you for improving this repository.

简述：欢迎补充文档、修复工作流问题、改进配置接口，但请不要把私有流域信息重新带回公开仓库。

## Ground Rules

- Keep the repository basin-agnostic.
- Do not add private raw data, private result snapshots, or case-study text.
- Keep public interfaces in English.
- Short Chinese notes are acceptable in documentation where they improve clarity for the intended audience.
- Preserve the output contracts used by the three model tiers unless there is a strong migration reason.

## Before You Open a Pull Request

1. Make sure new docs or code do not contain local drive-letter paths.
2. Make sure no case-study names, run IDs, or private filenames appear in committed content.
3. Keep configuration changes schema-compatible or document the migration.
4. Update relevant docs when you change workflow behavior.

## Documentation Changes

When changing docs:

- prefer concise, operational language
- state prerequisites and limitations clearly
- keep English primary, with concise Chinese support only where useful

## Code Changes

When changing workflow code:

- prefer configuration over hidden heuristics
- avoid hardcoded absolute paths
- keep optional dependencies optional
- document new CLI arguments in the README or docs

## Issue Reporting

Good issues include:

- the entrypoint you ran
- the config section involved
- the expected output
- the actual output or error
- whether the issue reproduces on the synthetic example

## Pull Request Checklist

- [ ] No private data or private run outputs added
- [ ] No basin-specific content added
- [ ] Docs updated if interface or behavior changed
- [ ] Output contracts preserved or explicitly documented
- [ ] New dependencies documented
