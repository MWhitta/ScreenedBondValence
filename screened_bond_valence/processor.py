"""File-oriented batch processor for ScreenedBondValence."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .api import ScreenedBondValenceService, fetch_materials_project_inputs
from .constants import DEFAULT_ALGORITHMS


class BondValenceProcessor:
    """Batch processor that writes the on-disk `res/` tree."""

    def __init__(
        self,
        api_key: str,
        algos: Sequence[str] | None,
        cations: Sequence[str],
        anions: Sequence[str],
        *,
        output_dir: str | Path = "res",
        r0_bounds: tuple[float, float] = (0.0, 5.0),
        energy_above_hull: tuple[float, float] = (0.0, 0.05),
    ):
        self.api_key = api_key
        self.algos = tuple(algos or DEFAULT_ALGORITHMS)
        self.cations = tuple(cations)
        self.anions = tuple(anions)
        self.output_dir = Path(output_dir)
        self.energy_above_hull = energy_above_hull
        self.service = ScreenedBondValenceService(algorithms=self.algos, r0_bounds=r0_bounds)
        self._ensure_directories()

    def _ensure_directories(self) -> None:
        self.output_dir.mkdir(exist_ok=True)
        for cation in self.cations:
            for anion in self.anions:
                pair_dir = self.output_dir / f"{cation}{anion}"
                (pair_dir / "params").mkdir(parents=True, exist_ok=True)
                (pair_dir / "R0Bs").mkdir(exist_ok=True)
                (pair_dir / "no_solu").mkdir(exist_ok=True)
                for algo in self.algos:
                    (pair_dir / "R0Bs" / algo).mkdir(exist_ok=True)

    def process_cation_system(self, cation: str, anion: str):
        print(f"Processing {cation}-{anion} system...")
        pair_dir = self.output_dir / f"{cation}{anion}"
        materials = fetch_materials_project_inputs(
            self.api_key,
            cation=cation,
            anion=anion,
            energy_above_hull=self.energy_above_hull,
        )
        self._save_possible_species(pair_dir, materials)

        if not materials:
            print(f"No materials with possible species found for {cation}-{anion}; skipping.")
            return []

        previous_results = self._load_previous_results(pair_dir)
        pending_materials = [
            material
            for material in materials
            if material.material_id not in previous_results["solved"]
        ]
        if not pending_materials:
            print(f"All {cation}-{anion} materials already processed; nothing to do.")
            return []

        if previous_results["solved"]:
            print(
                f"Resuming {cation}-{anion}: "
                f"skipping {len(previous_results['solved'])} previously solved materials."
            )

        results = self.service.fit_many(
            pending_materials,
            progress=True,
            desc=f"Processing {cation}-{anion} materials",
        )
        self._save_results(
            pair_dir,
            results,
            existing_sijs=previous_results["sij"],
            existing_charges=previous_results["charges"],
            existing_no_solution=previous_results["no_solution"],
        )
        return results

    def _save_possible_species(self, pair_dir: Path, materials) -> None:
        species_data = {
            material.material_id: list(material.possible_species)
            for material in materials
            if material.possible_species
        }
        output_file = pair_dir / "params" / "dict_matID_possible_species.json"
        output_file.write_text(json.dumps(species_data, indent=2), encoding="utf-8")

    def _load_previous_results(self, pair_dir: Path) -> dict[str, Any]:
        solved_sets: list[set[str]] = []
        for algorithm in self.algos:
            algo_dir = pair_dir / "R0Bs" / algorithm
            solved_files = {path.stem for path in algo_dir.glob("*.txt")}
            if solved_files:
                solved_sets.append(solved_files)

        no_solution_rows: list[tuple[str, str, str, str | None, str]] = []
        for algorithm in self.algos:
            no_solution_rows.extend(self._read_no_solution_rows(pair_dir / "no_solu" / f"{algorithm}.txt"))

        return {
            "solved": set.intersection(*solved_sets) if solved_sets else set(),
            "sij": self._read_json_dict(pair_dir / "dict_sijs.json"),
            "charges": self._read_json_dict(pair_dir / "dict_charges.json"),
            "no_solution": self._dedupe_no_solution_rows(no_solution_rows),
        }

    def _read_json_dict(self, path: Path) -> dict[str, Any]:
        if not path.exists():
            return {}
        return json.loads(path.read_text(encoding="utf-8"))

    def _read_no_solution_rows(self, path: Path) -> list[tuple[str, str, str, str | None, str]]:
        if not path.exists():
            return []

        rows: list[tuple[str, str, str, str | None, str]] = []
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.rstrip("\n")
                if not stripped:
                    continue
                parts = stripped.split("\t")
                parts += [""] * (5 - len(parts))
                material_id, cation, anion, formula_pretty, reason = parts[:5]
                rows.append(
                    (
                        material_id,
                        cation,
                        anion,
                        formula_pretty or None,
                        reason,
                    )
                )
        return rows

    def _dedupe_no_solution_rows(
        self,
        rows: Sequence[tuple[str, str, str, str | None, str]],
    ) -> list[tuple[str, str, str, str | None, str]]:
        deduped: list[tuple[str, str, str, str | None, str]] = []
        seen: set[tuple[str, str, str, str | None, str]] = set()
        for row in rows:
            if row in seen:
                continue
            seen.add(row)
            deduped.append(row)
        return deduped

    def _save_results(
        self,
        pair_dir: Path,
        results,
        *,
        existing_sijs: dict[str, dict[str, float]] | None = None,
        existing_charges: dict[str, dict[str, float]] | None = None,
        existing_no_solution: Sequence[tuple[str, str, str, str | None, str]] | None = None,
    ) -> None:
        dict_sijs: dict[str, dict[str, float]] = dict(existing_sijs or {})
        dict_charges: dict[str, dict[str, float]] = dict(existing_charges or {})
        no_solution_rows: list[tuple[str, str, str, str | None, str]] = list(existing_no_solution or [])

        for result in results:
            if result.theoretical is not None and result.theoretical.has_solution:
                dict_sijs[result.material.material_id] = dict(result.theoretical.bond_valences)
                dict_charges[result.material.material_id] = dict(result.theoretical.charge_map)

            if result.failure_reasons:
                for reason in result.failure_reasons:
                    no_solution_rows.append(
                        (
                            result.material.material_id,
                            result.material.cation,
                            result.material.anion,
                            result.material.formula_pretty,
                            reason,
                        )
                    )

            for fit in result.fits:
                output_file = pair_dir / "R0Bs" / fit.algorithm / f"{result.material.material_id}.txt"
                np.savetxt(output_file, np.asarray([fit.r0, fit.b]))

        (pair_dir / "dict_sijs.json").write_text(json.dumps(dict_sijs, indent=2), encoding="utf-8")
        (pair_dir / "dict_charges.json").write_text(json.dumps(dict_charges, indent=2), encoding="utf-8")
        no_solution_rows = self._dedupe_no_solution_rows(no_solution_rows)

        for algorithm in self.algos:
            output_file = pair_dir / "no_solu" / f"{algorithm}.txt"
            with output_file.open("w", encoding="utf-8") as handle:
                for row in no_solution_rows:
                    handle.write("\t".join("" if value is None else str(value) for value in row) + "\n")
