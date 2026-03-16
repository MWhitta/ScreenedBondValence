# ADR-0005: Keep Package And CLI Imports Lazy

## Status

Accepted

## Context

Some workflows only need parser or data-model access, while others require
heavier scientific or web dependencies. Eager imports make simple commands
slower and can fail in partially provisioned environments.

## Decision

Keep package exports and CLI web/batch loading lazy where practical.
Optional or heavy modules should only be imported when their workflow is
actually invoked.

## Consequences

- `python -m screened_bond_valence --help` stays lightweight.
- Batch/API workflows remain usable even if web-only dependencies are missing.
- New modules should avoid adding unnecessary top-level imports to shared
  entrypoints.
