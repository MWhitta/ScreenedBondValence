# ADR-0004: Use Compact Summary Payloads For Downstream Integrations

## Status

Accepted

## Context

Downstream repositories do not need all internal fitting artifacts. They need
stable, compact summaries that are easy to serialize and version.

## Decision

Downstream integrations should consume grouped outputs produced by
`build_summary_payload()` using a serializer function. The preferred stable
record shape is `MaterialFitSummary.to_compact_dict()`.

## Consequences

- Downstream repositories are insulated from internal dataclass changes.
- Integration contracts can be documented as JSON schemas instead of Python
  object graphs.
- If the compact record changes, the integration spec and golden tests must be
  updated together.
