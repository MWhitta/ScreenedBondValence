"""High-level reusable APIs for fitting screened bond valence parameters."""

from __future__ import annotations

import json
from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from pathlib import Path
from statistics import median
from typing import Any

from mp_api.client import MPRester
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core.structure import Structure
from tqdm import tqdm

from .models import (
    AlgorithmFitResult,
    BondValenceMaterial,
    MaterialFitResult,
    MaterialFitSummary,
    TheoreticalBondValenceResult,
)
from .constants import DEFAULT_ALGORITHMS
from .solvers import BVParamSolver, TheoreticalBondValenceSolver

SummarySerializer = Callable[[MaterialFitSummary], dict[str, Any]]


def _extract_charge_map_from_structure(structure: Structure) -> dict[str, float] | None:
    charge_map: dict[str, float] = {}
    for site in structure:
        oxi_state = getattr(site.specie, "oxi_state", None)
        if oxi_state is None:
            return None
        charge_map[site.specie.symbol] = float(oxi_state)
    return charge_map or None


def _normalize_charge_map(charge_map: Mapping[str, float] | None) -> dict[str, float] | None:
    if charge_map is None:
        return None
    return {str(element): float(charge) for element, charge in charge_map.items()}


def _guess_oxidation_state_decorated_structure(structure: Structure) -> tuple[Structure, str]:
    try:
        decorated = BVAnalyzer().get_oxi_state_decorated_structure(structure.copy())
        if _extract_charge_map_from_structure(decorated) is not None:
            return decorated, "bv_analyzer"
    except Exception:
        pass

    decorated = structure.copy()
    try:
        decorated.add_oxidation_state_by_guess()
    except Exception as exc:
        raise ValueError("Unable to infer oxidation states for CIF structure") from exc

    if _extract_charge_map_from_structure(decorated) is None:
        raise ValueError("Unable to infer oxidation states for CIF structure")
    return decorated, "composition_guess"


def _prepare_structure_with_oxidation_states(
    structure: Structure,
    charge_map: Mapping[str, float] | None = None,
) -> tuple[Structure, dict[str, float], str]:
    normalized_charge_map = _normalize_charge_map(charge_map)
    required_elements = {element.symbol for element in structure.composition.elements}

    if normalized_charge_map and required_elements.issubset(normalized_charge_map):
        decorated = structure.copy()
        decorated.add_oxidation_state_by_element(normalized_charge_map)
        return decorated, normalized_charge_map, "charge_map"

    existing_charge_map = _extract_charge_map_from_structure(structure)
    if existing_charge_map is not None:
        if normalized_charge_map:
            merged_charge_map = dict(existing_charge_map)
            merged_charge_map.update(normalized_charge_map)
            decorated = structure.copy()
            decorated.add_oxidation_state_by_element(merged_charge_map)
            return decorated, merged_charge_map, "charge_map_override"
        return structure, existing_charge_map, "cif"

    decorated, source = _guess_oxidation_state_decorated_structure(structure)
    guessed_charge_map = _extract_charge_map_from_structure(decorated)
    if guessed_charge_map is None:
        raise ValueError("Unable to infer oxidation states for CIF structure")

    if normalized_charge_map:
        guessed_charge_map.update(normalized_charge_map)
        decorated = structure.copy()
        decorated.add_oxidation_state_by_element(guessed_charge_map)
        return decorated, guessed_charge_map, "charge_map_override"

    return decorated, guessed_charge_map, source


def build_material_from_cif(
    cif_path: str | Path,
    *,
    cation: str,
    anion: str,
    material_id: str | None = None,
    possible_species: Sequence[str] = (),
    charge_map: Mapping[str, float] | None = None,
    crystal_nn_kwargs: Mapping[str, Any] | None = None,
    formula_pretty: str | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> BondValenceMaterial:
    path = Path(cif_path)
    structure = Structure.from_file(path)
    prepared_structure, extracted_charge_map, oxidation_state_source = _prepare_structure_with_oxidation_states(
        structure,
        charge_map=charge_map,
    )
    crystal_nn = CrystalNN(cation_anion=True, distance_cutoffs=[0, 1], **dict(crystal_nn_kwargs or {}))
    structure_graph = crystal_nn.get_bonded_structure(prepared_structure)

    resolved_metadata = {
        "source": "cif",
        "path": str(path),
        "oxidation_state_source": oxidation_state_source,
    }
    if metadata:
        resolved_metadata.update(dict(metadata))

    return BondValenceMaterial(
        material_id=material_id or path.stem,
        cation=cation,
        anion=anion,
        structure=prepared_structure,
        structure_graph=structure_graph,
        possible_species=possible_species,
        formula_pretty=formula_pretty or prepared_structure.composition.reduced_formula,
        charge_map=extracted_charge_map,
        metadata=resolved_metadata,
    )


def fetch_materials_project_inputs(
    api_key: str,
    *,
    cation: str,
    anion: str,
    energy_above_hull: tuple[float, float] = (0.0, 0.05),
    require_possible_species: bool = True,
) -> list[BondValenceMaterial]:
    with MPRester(api_key=api_key) as mpr:
        summary_docs = mpr.materials.summary.search(
            elements=[cation, anion],
            energy_above_hull=energy_above_hull,
            fields=["material_id", "possible_species"],
        )

    species_by_material: dict[str, tuple[str, ...]] = {}
    material_ids: list[str] = []
    for doc in summary_docs:
        material_id = str(doc.material_id)
        possible_species = tuple(doc.possible_species or ())
        if require_possible_species and not possible_species:
            continue
        species_by_material[material_id] = possible_species
        material_ids.append(material_id)

    if not material_ids:
        return []

    with MPRester(api_key=api_key) as mpr:
        bond_docs = mpr.materials.bonds.search(
            material_ids=material_ids,
            fields=["material_id", "structure_graph", "formula_pretty"],
        )

    materials: list[BondValenceMaterial] = []
    for doc in bond_docs:
        material_id = str(doc.material_id)
        materials.append(
            BondValenceMaterial(
                material_id=material_id,
                cation=cation,
                anion=anion,
                structure=doc.structure_graph.structure,
                structure_graph=doc.structure_graph,
                possible_species=species_by_material.get(material_id, ()),
                formula_pretty=getattr(doc, "formula_pretty", None),
                metadata={
                    "source": "materials_project",
                    "energy_above_hull": list(energy_above_hull),
                },
            )
        )

    return materials


class ScreenedBondValenceService:
    """Reusable in-memory API for screened bond valence fitting."""

    def __init__(
        self,
        *,
        algorithms: Sequence[str] = DEFAULT_ALGORITHMS,
        r0_bounds: tuple[float, float] = (0.0, 5.0),
        element2charge: str | Path | Mapping[str, float] | None = None,
    ):
        self.algorithms = tuple(algorithms)
        self.r0_bounds = r0_bounds
        self.theoretical_solver = TheoreticalBondValenceSolver(element2charge=element2charge)

    def compute_theoretical(self, material: BondValenceMaterial) -> TheoreticalBondValenceResult:
        structure = getattr(material.structure_graph, "structure", material.structure)
        bond_valences, bond_types, bond_lengths, charge_map = self.theoretical_solver.get_sij(
            material.material_id,
            structure,
            material.structure_graph,
            possible_species=material.possible_species,
            charge_map=material.charge_map,
        )
        return TheoreticalBondValenceResult(
            material_id=material.material_id,
            bond_valences=bond_valences or {},
            bond_types=bond_types,
            bond_lengths=bond_lengths,
            charge_map=charge_map,
        )

    def fit_material(self, material: BondValenceMaterial) -> MaterialFitResult:
        theoretical = self.compute_theoretical(material)
        if not theoretical.has_solution:
            return MaterialFitResult(
                material=material,
                theoretical=theoretical,
                failure_reasons=("no_network_solution",),
            )

        fits: list[AlgorithmFitResult] = []
        failure_reasons: list[str] = []

        for algorithm in self.algorithms:
            solver = BVParamSolver(algo=algorithm)
            solution = solver.solve_R0Bs(
                cation=material.cation,
                anion=material.anion,
                bond_type_list=theoretical.bond_types,
                networkValence_dict=theoretical.bond_valences,
                bondLen_dict=theoretical.bond_lengths,
                materID=material.material_id,
                chem_formula=material.formula_pretty,
                R0_bounds=self.r0_bounds,
            )
            if solution is None:
                if solver.last_failure_reason:
                    failure_reasons.append(f"{algorithm}:{solver.last_failure_reason}")
                continue

            fits.append(AlgorithmFitResult(algorithm=algorithm, r0=solution[0], b=solution[1]))

        return MaterialFitResult(
            material=material,
            theoretical=theoretical,
            fits=fits,
            failure_reasons=tuple(dict.fromkeys(failure_reasons)),
        )

    def fit_many(
        self,
        materials: Iterable[BondValenceMaterial],
        *,
        progress: bool = False,
        desc: str = "Fitting bond valence parameters",
    ) -> list[MaterialFitResult]:
        material_list = list(materials)
        iterator = material_list
        if progress:
            iterator = tqdm(material_list, desc=desc)
        return [self.fit_material(material) for material in iterator]

    def summarize(self, result: MaterialFitResult, reducer: Callable[[Sequence[float]], float] = median) -> MaterialFitSummary | None:
        return result.aggregate(reducer=reducer)


def build_summary_payload(
    results: Iterable[MaterialFitResult],
    *,
    classifier: Callable[[MaterialFitSummary], str | None] | None = None,
    reducer: Callable[[Sequence[float]], float] = median,
    serializer: SummarySerializer | None = None,
) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    key_name = "results"
    serialize = serializer or (lambda summary: summary.to_dict())

    for result in results:
        summary = result.aggregate(reducer=reducer)
        if summary is None:
            continue

        key = key_name if classifier is None else classifier(summary)
        if key is None:
            continue
        grouped[key].append(serialize(summary))

    return dict(grouped)


def export_summary_payload(
    output_path: str | Path,
    results: Iterable[MaterialFitResult],
    *,
    classifier: Callable[[MaterialFitSummary], str | None] | None = None,
    reducer: Callable[[Sequence[float]], float] = median,
    serializer: SummarySerializer | None = None,
) -> dict[str, list[dict[str, Any]]]:
    payload = build_summary_payload(
        results,
        classifier=classifier,
        reducer=reducer,
        serializer=serializer,
    )
    Path(output_path).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload
