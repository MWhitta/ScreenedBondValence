"""Shared data models for the ScreenedBondValence API."""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import median, pstdev
from typing import Any, Callable, Mapping, Sequence


Reducer = Callable[[Sequence[float]], float]


def _copy_metadata(metadata: Mapping[str, Any] | None) -> dict[str, Any]:
    return dict(metadata or {})


@dataclass(slots=True)
class BondValenceMaterial:
    """All inputs required to fit one cation-anion system in one material."""

    material_id: str
    cation: str
    anion: str
    structure: Any
    structure_graph: Any
    possible_species: Sequence[str] = ()
    formula_pretty: str | None = None
    charge_map: Mapping[str, float] | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.possible_species = tuple(self.possible_species)
        self.metadata = _copy_metadata(self.metadata)
        if self.charge_map is not None:
            self.charge_map = {
                str(element): float(charge)
                for element, charge in self.charge_map.items()
            }


@dataclass(slots=True)
class TheoreticalBondValenceResult:
    """The solved network bond valences for a material."""

    material_id: str
    bond_valences: dict[str, float]
    bond_types: Sequence[str]
    bond_lengths: dict[str, float]
    charge_map: dict[str, float]

    def __post_init__(self) -> None:
        self.bond_types = tuple(self.bond_types)

    @property
    def has_solution(self) -> bool:
        return bool(self.bond_valences)

    def to_dict(self) -> dict[str, Any]:
        return {
            "material_id": self.material_id,
            "bond_valences": dict(self.bond_valences),
            "bond_types": list(self.bond_types),
            "bond_lengths": dict(self.bond_lengths),
            "charge_map": dict(self.charge_map),
        }


@dataclass(slots=True)
class AlgorithmFitResult:
    """One optimizer's fitted bond-valence parameters."""

    algorithm: str
    r0: float
    b: float

    def to_dict(self) -> dict[str, Any]:
        return {"algorithm": self.algorithm, "R0": self.r0, "B": self.b}


@dataclass(slots=True)
class MaterialFitSummary:
    """Aggregated per-material summary across one or more optimizers."""

    material_id: str
    cation: str
    anion: str
    formula_pretty: str | None
    r0: float
    b: float
    r0_std: float
    b_std: float
    n_algos: int
    algorithms: Sequence[str]
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.algorithms = tuple(self.algorithms)
        self.metadata = _copy_metadata(self.metadata)

    def to_dict(self) -> dict[str, Any]:
        payload = {
            "material_id": self.material_id,
            "mid": self.material_id,
            "cation": self.cation,
            "anion": self.anion,
            "formula_pretty": self.formula_pretty,
            "formula": self.formula_pretty,
            "R0": self.r0,
            "B": self.b,
            "R0_std": self.r0_std,
            "B_std": self.b_std,
            "n_algos": self.n_algos,
            "algorithms": list(self.algorithms),
        }
        if self.metadata:
            payload["metadata"] = dict(self.metadata)
        return payload

    def to_compact_dict(self) -> dict[str, Any]:
        return {
            "mid": self.material_id,
            "formula": self.formula_pretty,
            "R0": self.r0,
            "B": self.b,
            "R0_std": self.r0_std,
            "B_std": self.b_std,
            "n_algos": self.n_algos,
        }


@dataclass(slots=True)
class MaterialFitResult:
    """All fitting artifacts for one material."""

    material: BondValenceMaterial
    theoretical: TheoreticalBondValenceResult | None
    fits: Sequence[AlgorithmFitResult] = ()
    failure_reasons: Sequence[str] = ()

    def __post_init__(self) -> None:
        self.fits = tuple(self.fits)
        self.failure_reasons = tuple(self.failure_reasons)

    @property
    def has_fit(self) -> bool:
        return bool(self.fits)

    def aggregate(self, reducer: Reducer = median) -> MaterialFitSummary | None:
        if not self.fits:
            return None

        r0_values = [fit.r0 for fit in self.fits]
        b_values = [fit.b for fit in self.fits]
        return MaterialFitSummary(
            material_id=self.material.material_id,
            cation=self.material.cation,
            anion=self.material.anion,
            formula_pretty=self.material.formula_pretty,
            r0=float(reducer(r0_values)),
            b=float(reducer(b_values)),
            r0_std=float(pstdev(r0_values)),
            b_std=float(pstdev(b_values)),
            n_algos=len(self.fits),
            algorithms=[fit.algorithm for fit in self.fits],
            metadata=self.material.metadata,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "material_id": self.material.material_id,
            "cation": self.material.cation,
            "anion": self.material.anion,
            "formula_pretty": self.material.formula_pretty,
            "possible_species": list(self.material.possible_species),
            "metadata": dict(self.material.metadata),
            "charge_map": None if self.material.charge_map is None else dict(self.material.charge_map),
            "theoretical": None if self.theoretical is None else self.theoretical.to_dict(),
            "fits": [fit.to_dict() for fit in self.fits],
            "failure_reasons": list(self.failure_reasons),
        }
