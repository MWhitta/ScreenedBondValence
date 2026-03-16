# ADR-0003: Preserve A Resumable On-Disk Batch Output Tree

## Status

Accepted

## Context

Materials Project batch runs can be large and expensive. Recomputing all
materials from scratch on every interruption is wasteful.

## Decision

The batch processor will continue to write and resume from a `res/` tree with
per-algorithm fit files plus shared `dict_sijs.json`, `dict_charges.json`, and
`no_solu/` records.

Materials are considered solved only when outputs exist for every configured
algorithm.

## Consequences

- Long runs can resume without refitting completed materials.
- The on-disk layout remains inspectable and script-friendly.
- Any future output-layout change must preserve resumability semantics or
  explicitly replace them.
