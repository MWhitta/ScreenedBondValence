# ADR-0002: Resolve CIF Oxidation States With A Layered Policy

## Status

Accepted

## Context

Many CIF files do not carry complete oxidation state information, but the
bond-valence workflow requires a charge map. Blindly falling back to a static
element table hides structure-specific chemistry.

## Decision

For CIF workflows, oxidation states are resolved in this order:

1. explicit oxidation states already present in the CIF
2. `pymatgen.analysis.bond_valence.BVAnalyzer`
3. `Structure.add_oxidation_state_by_guess()`
4. explicit caller-provided charge overrides

The chosen source is recorded in material metadata.

## Consequences

- CIF fitting remains usable when oxidation states are omitted.
- Structure-aware guessing is preferred over generic defaults.
- Downstream consumers can inspect metadata to see how charges were assigned.
