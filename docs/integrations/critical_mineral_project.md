# Integration Spec: CriticalMineralProject

## Purpose

This document defines the current integration contract between
`ScreenedBondValence` and `CriticalMineralProject`.

## Producer

`ScreenedBondValence` produces grouped summary payloads via:

- `build_summary_payload(results, classifier=..., serializer=...)`
- `export_summary_payload(output_path, results, classifier=..., serializer=...)`

For `CriticalMineralProject`, the serializer should be:

```python
serializer=lambda summary: summary.to_compact_dict()
```

## Consumer Expectations

`CriticalMineralProject` currently expects JSON shaped like:

```json
{
  "oxides": [
    {
      "mid": "mp-1019547",
      "formula": "BaLi2Mg(PO4)2",
      "R0": 0.0,
      "B": 1.5125929975078756,
      "R0_std": 1.1847948628290317,
      "B_std": 0.8546488334257794,
      "n_algos": 5
    }
  ],
  "hydroxides": []
}
```

## Required Top-Level Shape

- JSON object
- keys are downstream-defined buckets such as `oxides` and `hydroxides`
- values are arrays of compact summary records

## Required Record Fields

- `mid`
  Materials Project ID or other material identifier
- `formula`
  reduced or pretty formula used for labeling
- `R0`
  aggregated bond-valence `R0`
- `B`
  aggregated bond-valence `B`
- `R0_std`
  standard deviation across optimizer outputs
- `B_std`
  standard deviation across optimizer outputs
- `n_algos`
  number of optimizer results contributing to the aggregate

## Current Aggregation Policy

- per-material aggregation uses `MaterialFitResult.aggregate()`
- the default reducer is the median across optimizer outputs
- downstream grouping logic is owned by `CriticalMineralProject`, not by
  `ScreenedBondValence`

## Compatibility Rules

- Adding new top-level buckets is allowed if the consumer expects them.
- Adding extra per-record keys is allowed only if existing required keys remain
  unchanged.
- Renaming or removing the required fields above is a breaking change.
- Any intentional contract change must update this spec, `README.md`, and the
  golden tests.
