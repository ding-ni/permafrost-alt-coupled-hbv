# Architecture

This repository is structured as a public workflow layer around a minimal vendored legacy core.

## Main Layers

### 1. Public Entry Layer

Located in `workflow/`.

This is the main interface surface for users. It provides:

- preprocessing entrypoints
- model entrypoints
- analysis entrypoints

These scripts should be treated as the stable public surface of the repository.

### 2. Shared Runtime Layer

Located in `shared/`.

This layer centralizes:

- config loading
- workspace path construction
- public-to-legacy compatibility helpers
- parameter parsing and output handoff helpers

This keeps wrapper behavior consistent across the workflow.

### 3. Legacy Core Layer

Located in `legacy_core/`.

This is the vendored minimal execution core retained for reproducibility and interface continuity. It exists so the public wrappers do not depend on a private project workspace layout.

## Design Intent

The repository is intentionally organized so that:

- users interact with `workflow/`
- maintainers evolve public interfaces in `shared/`
- low-level retained execution code remains isolated inside `legacy_core/`

## Why This Matters

This separation makes it easier to:

- document the public interface clearly
- keep private historical structure out of the published repository
- eventually refactor or replace legacy internals without breaking the public surface all at once
