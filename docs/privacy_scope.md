# Privacy and Publication Scope

This repository is a public method release.

简述：公开版只保留可复用的方法，不保留任何特定流域的建模叙述、默认案例或结果快照。

## Excluded Material

The public release must not include:

- case-study names
- thesis chapter text
- private configuration examples
- private run directories
- local drive-letter paths
- figure captions or labels tied to a private study basin

## Allowed Material

The public release may include:

- method descriptions
- generic diagrams
- schema definitions
- workflow entrypoints
- synthetic example data meant only for smoke testing

## Maintainer Check

Before publishing, search the repository for:

- basin-specific names
- local absolute paths
- private run IDs
- private data filenames

Any match should be treated as a release blocker unless it is part of a clearly documented placeholder.
