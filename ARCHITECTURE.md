# Architecture

## Overview

`ScreenedBondValence` computes theoretical bond valences from crystal
structures and fits bond-valence parameters `R0` and `B` for cation-anion
pairs. The repository is organized as a single Python package with one CLI and
one offline test suite.

## Package Layout

- `screened_bond_valence/api.py`
  High-level in-memory workflows. This is the main entrypoint for building
  materials from CIFs, fetching Materials Project inputs, fitting materials,
  and exporting grouped payloads.
- `screened_bond_valence/solvers.py`
  Low-level network and parameter solvers. This module computes `Sij` values
  from valence-sum and loop equations, then fits `R0/B` using SciPy global
  optimizers.
- `screened_bond_valence/processor.py`
  File-oriented batch workflow. It fetches Materials Project data, resumes
  previous runs from the on-disk `res/` tree, and writes fit outputs.
- `screened_bond_valence/models.py`
  Shared data models for inputs, intermediate results, and aggregated fit
  summaries.
- `screened_bond_valence/web.py`
  Gradio UI for fitting one CIF interactively.
- `screened_bond_valence/cli.py`
  Command-line interface for batch and web workflows.
- `screened_bond_valence/constants.py`
  Lightweight shared constants such as the supported optimizer list.

## Main Data Flow

### CIF Workflow

1. `build_material_from_cif()` loads a structure from disk.
2. Oxidation states are resolved in priority order:
   explicit CIF values, `BVAnalyzer`, composition-based guess, then explicit
   overrides if provided.
3. `CrystalNN` builds the bonded structure graph.
4. `ScreenedBondValenceService.fit_material()` computes theoretical `Sij`
   values and runs the configured optimizers.
5. The result can be aggregated with `MaterialFitResult.aggregate()` and
   serialized with `build_summary_payload()`.

### Materials Project Batch Workflow

1. `fetch_materials_project_inputs()` retrieves summary and bonding docs.
2. `BondValenceProcessor.process_cation_system()` writes
   `dict_matID_possible_species.json`.
3. Existing `R0Bs`, `dict_sijs.json`, `dict_charges.json`, and `no_solu/`
   files are loaded to determine resumable state.
4. Only unsolved materials are fitted.
5. New outputs are merged back into the `res/<cation><anion>/` tree.

## Core Models

- `BondValenceMaterial`
  All information needed to fit one material for one cation-anion pair.
- `TheoreticalBondValenceResult`
  Solved network bond valences, bond lengths, and charge map.
- `MaterialFitResult`
  One material’s theoretical result plus all optimizer outputs.
- `MaterialFitSummary`
  Aggregated per-material `R0/B` summary used for downstream payloads.

## Outputs

### Batch Output Tree

The batch processor writes:

- `params/dict_matID_possible_species.json`
- `R0Bs/<algorithm>/<material_id>.txt`
- `no_solu/<algorithm>.txt`
- `dict_sijs.json`
- `dict_charges.json`

### Downstream JSON Payload

Downstream integrations should prefer `build_summary_payload()` with a custom
classifier and serializer. The compact serializer
`MaterialFitSummary.to_compact_dict()` is the current contract used by
`CriticalMineralProject`.

## Stability Boundaries

- `screened_bond_valence.api`, `screened_bond_valence.processor`, and the CLI
  are the supported workflow surfaces.
- `screened_bond_valence.solvers` is more implementation-oriented and may
  change as the fitting approach evolves.
- The compact summary payload is treated as a stable integration contract.

## Known Technical Debt

- The `R0/B` fitting stage still uses global optimizers for a two-parameter
  problem. A future least-squares reformulation is a likely improvement.
- Theoretical `Sij` solving depends on SymPy equation solving and can be slow
  for larger systems.
- The sample CIF coverage is narrow. Additional representative fixtures would
  improve future regression detection.
