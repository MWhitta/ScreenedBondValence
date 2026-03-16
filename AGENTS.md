# AGENTS.md

This file is for human contributors and LLM-based coding agents working in
`ScreenedBondValence`.

## Repository Rules

- Keep the repository package-first. New functionality belongs under
  `screened_bond_valence/`, not as top-level scripts.
- Treat `pyproject.toml` as the single source of truth for runtime
  dependencies. `requirements.txt` should remain a thin editable-install shim.
- Do not commit generated artifacts such as `__pycache__/`, notebook outputs,
  `res/`, or local virtual environments.
- Keep offline tests deterministic. Unit tests must not require network
  access or a Materials Project API key.
- Prefer stable public contracts over prompt history. If behavior changes,
  update the docs and tests rather than adding ad hoc context notes.

## Workflows

- Install: `pip install -e .`
- Run tests: `PYTHONPATH=. python -m unittest discover -s tests -v`
- Batch CLI: `screened-bond-valence batch ...`
- Web CLI: `screened-bond-valence web ...`

## Change Checklist

Update the following when relevant:

- `README.md` for user-facing workflow changes
- `ARCHITECTURE.md` for package structure or data-flow changes
- `docs/adr/*.md` when a design decision changes
- `docs/integrations/critical_mineral_project.md` when the downstream JSON
  contract changes
- `tests/golden/*` and `tests/test_golden_outputs.py` when stable outputs are
  intentionally changed

## Current Constraints

- Oxidation states for CIF workflows are resolved in this order:
  explicit CIF values, `BVAnalyzer`, then
  `Structure.add_oxidation_state_by_guess()`, with user overrides allowed.
- Batch processing must remain resumable using the existing `res/` tree.
- The compact summary payload is the preferred integration surface for
  downstream projects such as `CriticalMineralProject`.
- The CLI and package imports should stay lightweight. Avoid eager imports of
  optional or heavy modules when they are only needed for specific commands.

## Preferred Extension Points

- Add new high-level workflows in `api.py`
- Add or change fitting/network logic in `solvers.py`
- Extend on-disk batch orchestration in `processor.py`
- Keep web-only code in `web.py`
- Expose new public APIs through `screened_bond_valence/__init__.py`
