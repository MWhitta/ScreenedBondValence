"""Reusable APIs for ScreenedBondValence."""

from __future__ import annotations

from importlib import import_module


_EXPORTS = {
    "AlgorithmFitResult": ("screened_bond_valence.models", "AlgorithmFitResult"),
    "BondValenceMaterial": ("screened_bond_valence.models", "BondValenceMaterial"),
    "BondValenceProcessor": ("screened_bond_valence.processor", "BondValenceProcessor"),
    "BVParamSolver": ("screened_bond_valence.solvers", "BVParamSolver"),
    "DEFAULT_ALGORITHMS": ("screened_bond_valence.constants", "DEFAULT_ALGORITHMS"),
    "MaterialFitResult": ("screened_bond_valence.models", "MaterialFitResult"),
    "MaterialFitSummary": ("screened_bond_valence.models", "MaterialFitSummary"),
    "ScreenedBondValenceService": ("screened_bond_valence.api", "ScreenedBondValenceService"),
    "TheoreticalBondValenceResult": ("screened_bond_valence.models", "TheoreticalBondValenceResult"),
    "TheoreticalBondValenceSolver": ("screened_bond_valence.solvers", "TheoreticalBondValenceSolver"),
    "build_material_from_cif": ("screened_bond_valence.api", "build_material_from_cif"),
    "build_summary_payload": ("screened_bond_valence.api", "build_summary_payload"),
    "export_summary_payload": ("screened_bond_valence.api", "export_summary_payload"),
    "fetch_materials_project_inputs": ("screened_bond_valence.api", "fetch_materials_project_inputs"),
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    try:
        module_name, attribute_name = _EXPORTS[name]
    except KeyError as exc:  # pragma: no cover - Python import protocol
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc

    value = getattr(import_module(module_name), attribute_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
