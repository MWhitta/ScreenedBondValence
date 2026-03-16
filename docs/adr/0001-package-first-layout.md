# ADR-0001: Use A Package-First Layout With One CLI

## Status

Accepted

## Context

The repository originally exposed multiple top-level scripts with overlapping
responsibilities. That made packaging, testing, and future extension harder.

## Decision

The repository will keep one installable package under
`screened_bond_valence/` and one CLI entrypoint,
`screened-bond-valence`.

## Consequences

- Public workflows live in importable modules rather than top-level scripts.
- Tests can target importable functions and classes directly.
- Packaging and downstream reuse are simpler.
- New features should not be added as standalone root-level scripts.
