"""Low-level network and parameter solvers for ScreenedBondValence."""

from __future__ import annotations

import json
import os
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from scipy.optimize import brute, differential_evolution, direct, dual_annealing, shgo


DEFAULT_ELEMENT2CHARGE_PATH = Path(__file__).resolve().parents[1] / "element2charge.json"


def _load_mapping(
    source: str | Path | Mapping[str, Any] | None,
    *,
    default: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if source is None:
        return dict(default or {})
    if isinstance(source, Mapping):
        return dict(source)

    with Path(source).open("r", encoding="utf-8") as handle:
        return json.load(handle)


class TheoreticalBondValenceSolver:
    """Solve theoretical bond valences from valence-sum and loop equations."""

    def __init__(
        self,
        element2charge: str | Path | Mapping[str, float] | None = None,
        species_by_material: str | Path | Mapping[str, Sequence[str]] | None = None,
    ):
        self.dict_ele2charge = _load_mapping(
            element2charge or DEFAULT_ELEMENT2CHARGE_PATH,
        )
        self.dict_species_matID = _load_mapping(species_by_material, default={})

    def get_element(self, label: str) -> str:
        return re.split(r"[^a-zA-Z]", label)[0]

    def idx2label(self, struct: Structure) -> dict[int, str]:
        return {i: site.label for i, site in enumerate(struct.sites)}

    def relabel_site_labels_sym(self, struct: Structure) -> None:
        sga = SpacegroupAnalyzer(struct)
        sym_struct = sga.get_symmetrized_structure()

        new_labels: dict[int, str] = {}
        processed_element: dict[str, int] = {}

        for group in sym_struct.equivalent_indices:
            species = self.get_element(struct[group[0]].species_string)
            processed_element[species] = processed_element.get(species, 0) + 1
            for site in group:
                new_labels[site] = f"{species}{processed_element[species]}"

        for i, site in enumerate(struct):
            site.label = new_labels[i]

    def get_element2charge(self, list_of_possible_species: Sequence[str]) -> dict[str, float]:
        ele2char: dict[str, float] = {}
        for spec in list_of_possible_species:
            element = re.split(r"[^a-zA-Z]", spec)[0]
            charge = spec.split(element)[1]
            chargeval = re.split(r"[+,-]", charge)[0] or "1"
            sign = charge[-1]
            ele2char[element] = float(f"{sign}{chargeval}")
        return ele2char

    def resolve_charge_map(
        self,
        mat_id: str,
        *,
        possible_species: Sequence[str] | None = None,
        charge_map: Mapping[str, float] | None = None,
    ) -> dict[str, float]:
        if charge_map:
            return {str(element): float(charge) for element, charge in charge_map.items()}

        species = tuple(possible_species or self.dict_species_matID.get(mat_id, ()))
        return self.get_element2charge(species) or dict(self.dict_ele2charge)

    def find_cycles(self, edges: Sequence[tuple[int, int]]) -> list[tuple[int, ...]]:
        graph: dict[int, set[int]] = {}
        for node1, node2 in edges:
            graph.setdefault(node1, set()).add(node2)
            graph.setdefault(node2, set()).add(node1)

        cycles: set[tuple[int, ...]] = set()

        def dfs(node: int, visited: set[int], path: list[int]) -> None:
            visited.add(node)
            path.append(node)

            for neighbor in graph[node]:
                if neighbor in path:
                    idx = path.index(neighbor)
                    cycles.add(tuple(path[idx:]))
                elif neighbor not in visited:
                    dfs(neighbor, visited, path)

            path.pop()

        visited: set[int] = set()
        for node in graph:
            if node not in visited:
                dfs(node, visited, [])

        return [cycle for cycle in cycles if len(cycle) >= 4]

    def get_cycle_path(self, graph: Any, idx2label: Mapping[int, str]) -> list[list[str]]:
        cycles = self.find_cycles(list(graph.edges()))
        return [[idx2label[idx] for idx in cycle] for cycle in cycles]

    def get_eq_cycle(self, cur_cycle: Sequence[str], dict_element_charge: Mapping[str, float]) -> str:
        eq_temp = ""
        for i, label in enumerate(cur_cycle):
            atom_element = self.get_element(label)
            next_label = cur_cycle[0] if i == len(cur_cycle) - 1 else cur_cycle[i + 1]

            if dict_element_charge[atom_element] > 0:
                eq_temp = eq_temp[:-1]
                temp_bond = f"-{label}{next_label}+"
            else:
                temp_bond = f"{next_label}{label}+"
            eq_temp += temp_bond

        return eq_temp[:-1]

    def solve_sij(
        self,
        variables: Sequence[str],
        equations: Sequence[tuple[str, float]],
    ) -> dict[str, float] | None:
        symbols_list = sp.symbols(list(variables))
        eqs = [sp.Eq(sp.sympify(eq[0]), eq[1]) for eq in equations]
        solution = sp.solve(eqs, symbols_list)

        if not solution or any(isinstance(value, sp.Add) for value in solution.values()):
            return None

        return {str(symbol): float(value) for symbol, value in solution.items()}

    def get_eqs_from_valence_sum_rule(
        self,
        struct: Structure,
        bonds: Any,
        dict_charge: Mapping[str, float],
    ) -> tuple[list[tuple[str, float]], set[str], dict[str, float]]:
        equations: list[tuple[str, float]] = []
        bond_variables: set[str] = set()
        bond_lengths: dict[str, float] = {}

        for i, site in enumerate(struct):
            site_element = self.get_element(site.label)
            atom_valence = getattr(site.specie, "oxi_state", dict_charge[site_element])

            bond_terms: list[str] = []
            for neighbor in bonds.get_connected_sites(i):
                bond_type = (
                    f"{site.label}{neighbor.site.label}"
                    if atom_valence > 0
                    else f"{neighbor.site.label}{site.label}"
                )

                bond_variables.add(bond_type)
                bond_lengths[bond_type] = neighbor.dist
                bond_terms.append(bond_type)

            if bond_terms:
                equation = " + ".join(bond_terms)
                equations.append((equation, abs(atom_valence)))

        return equations, bond_variables, bond_lengths

    def get_eqs_from_loops(
        self,
        struct: Structure,
        graph: Any,
        dict_charge: Mapping[str, float],
    ) -> list[tuple[str, float]]:
        id2label = self.idx2label(struct)
        cycle_list = self.get_cycle_path(graph.graph, id2label)
        return [(self.get_eq_cycle(cur_cycle, dict_charge), 0) for cur_cycle in cycle_list]

    def get_sij(
        self,
        mat_id: str,
        struct: Structure,
        graph: Any,
        *,
        possible_species: Sequence[str] | None = None,
        charge_map: Mapping[str, float] | None = None,
    ) -> tuple[dict[str, float] | None, list[str], dict[str, float], dict[str, float]]:
        self.relabel_site_labels_sym(struct)
        dict_charge = self.resolve_charge_map(
            mat_id,
            possible_species=possible_species,
            charge_map=charge_map,
        )

        equations_val_sum, bond_vars, bond_lengths = self.get_eqs_from_valence_sum_rule(
            struct,
            graph,
            dict_charge,
        )
        equations_cycle = self.get_eqs_from_loops(struct, graph, dict_charge)

        sijs = self.solve_sij(bond_vars, equations_val_sum + equations_cycle)
        return sijs, list(bond_vars), bond_lengths, dict_charge


class BVParamSolver:
    """Optimize bond valence parameters for one material and bond family."""

    def __init__(self, save_dir: str | Path | None = None, algo: str = "shgo", no_sol: list | None = None):
        self.algo = algo
        self.no_sol = [] if no_sol is None else no_sol
        self.last_failure_reason: str | None = None

        if save_dir is not None:
            os.makedirs(Path(save_dir) / "R0Bs" / algo, exist_ok=True)
            os.makedirs(Path(save_dir) / "no_solu", exist_ok=True)

    def objective(self, variables: Sequence[float], eqs: Sequence[sp.Expr]) -> float:
        r0, b = sp.symbols("R0 B")
        return sum(float(eq.subs({r0: variables[0], b: variables[1]})) ** 2 for eq in eqs)

    def get_eqs_for_R0B(
        self,
        cation: str,
        anion: str,
        bond_type_list: Sequence[str],
        network_valence_dict: Mapping[str, float],
        bond_len_dict: Mapping[str, float],
        mat_id: str,
        reduced_formula: str | None,
        r0_bounds: tuple[float, float],
    ) -> tuple[list[str], tuple[float, float] | None]:
        self.last_failure_reason = None
        target_bonds = [
            bond_type
            for bond_type in bond_type_list
            if re.split(r"\d+", bond_type)[0] == cation
            and re.split(r"\d+", bond_type)[1] == anion
        ]

        bmax, bmin = -10.0, 1000.0
        eqs_list: list[str] = []

        for bond in target_bonds:
            sij = network_valence_dict[bond]
            if sij <= 0:
                self.last_failure_reason = "negative_Sij"
                self.no_sol.append((mat_id, cation, anion, reduced_formula, self.last_failure_reason))
                return [], None

            eqs_list.append(f"R0 - B*log({sij}) - {bond_len_dict[bond]}")

            if sij != 1:
                b1 = (r0_bounds[0] - bond_len_dict[bond]) / np.log(sij)
                b2 = (r0_bounds[1] - bond_len_dict[bond]) / np.log(sij)
                bmax = max(bmax, b1, b2)
                bmin = min(bmin, b1, b2)

        if not (-10 < bmax < np.inf) or not (-10 < bmin < np.inf):
            bmin, bmax = -5.0, 5.0

        if not eqs_list:
            self.last_failure_reason = "no_eqs_from_graph"
            self.no_sol.append((mat_id, cation, anion, reduced_formula, self.last_failure_reason))
            return [], None

        return eqs_list, (bmin, bmax)

    def solve_R0Bs(
        self,
        cation: str,
        anion: str,
        bond_type_list: Sequence[str],
        networkValence_dict: Mapping[str, float],
        bondLen_dict: Mapping[str, float],
        materID: str,
        chem_formula: str | None,
        R0_bounds: tuple[float, float],
    ) -> tuple[float, float] | None:
        self.last_failure_reason = None
        eqs_list_math, b_bounds = self.get_eqs_for_R0B(
            cation,
            anion,
            bond_type_list,
            networkValence_dict,
            bondLen_dict,
            materID,
            chem_formula,
            R0_bounds,
        )
        if not eqs_list_math or b_bounds is None:
            return None

        eqs_list_sympy = [sp.sympify(eq) for eq in eqs_list_math]
        optimizers = {
            "brute": lambda: brute(self.objective, [R0_bounds, b_bounds], args=(eqs_list_sympy,)),
            "shgo": lambda: shgo(self.objective, [R0_bounds, b_bounds], args=(eqs_list_sympy,)).x,
            "diff": lambda: differential_evolution(
                self.objective,
                [R0_bounds, b_bounds],
                args=(eqs_list_sympy,),
            ).x,
            "dual_annealing": lambda: dual_annealing(
                self.objective,
                [R0_bounds, b_bounds],
                args=(eqs_list_sympy,),
            ).x,
            "direct": lambda: direct(self.objective, [R0_bounds, b_bounds], args=(eqs_list_sympy,)).x,
        }

        algo = self.algo.lower()
        if algo not in optimizers:
            algo = "shgo"

        try:
            result = optimizers[algo]()
        except Exception as exc:  # pragma: no cover - optimizer-specific failure path
            self.last_failure_reason = f"optimizer_failed:{exc}"
            self.no_sol.append((materID, cation, anion, chem_formula, self.last_failure_reason))
            return None

        if isinstance(result, np.ndarray):
            return float(result[0]), float(result[1])
        if isinstance(result, (list, tuple)):
            return float(result[0]), float(result[1])
        return result

