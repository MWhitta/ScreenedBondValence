from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from screened_bond_valence.models import (
    AlgorithmFitResult,
    BondValenceMaterial,
    MaterialFitResult,
    TheoreticalBondValenceResult,
)
from screened_bond_valence.processor import BondValenceProcessor


def _make_material(material_id: str) -> BondValenceMaterial:
    return BondValenceMaterial(
        material_id=material_id,
        cation="Li",
        anion="O",
        structure=object(),
        structure_graph=object(),
        possible_species=("Li+", "O2-"),
        formula_pretty="Li2O",
    )


def _make_result(material_id: str) -> MaterialFitResult:
    material = _make_material(material_id)
    theoretical = TheoreticalBondValenceResult(
        material_id=material_id,
        bond_valences={"Li1O1": 1.0},
        bond_types=("Li1O1",),
        bond_lengths={"Li1O1": 2.0},
        charge_map={"Li": 1.0, "O": -2.0},
    )
    return MaterialFitResult(
        material=material,
        theoretical=theoretical,
        fits=(AlgorithmFitResult("shgo", 1.3, 0.5), AlgorithmFitResult("diff", 1.3, 0.5)),
    )


class BondValenceProcessorTests(unittest.TestCase):
    def test_resumes_and_merges_existing_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            processor = BondValenceProcessor(
                api_key="test-key",
                algos=("shgo", "diff"),
                cations=("Li",),
                anions=("O",),
                output_dir=tmpdir,
            )
            pair_dir = Path(tmpdir) / "LiO"
            for algorithm in ("shgo", "diff"):
                (pair_dir / "R0Bs" / algorithm / "mp-solved.txt").write_text("1.2\n0.4\n", encoding="utf-8")
                (pair_dir / "no_solu" / f"{algorithm}.txt").write_text(
                    "mp-failed\tLi\tO\tLiO\tno_network_solution\n",
                    encoding="utf-8",
                )

            (pair_dir / "dict_sijs.json").write_text(
                json.dumps({"mp-solved": {"Li1O1": 1.0}}),
                encoding="utf-8",
            )
            (pair_dir / "dict_charges.json").write_text(
                json.dumps({"mp-solved": {"Li": 1.0, "O": -2.0}}),
                encoding="utf-8",
            )

            materials = [_make_material("mp-solved"), _make_material("mp-new")]
            new_result = _make_result("mp-new")
            processor.service.fit_many = Mock(return_value=[new_result])

            with patch("screened_bond_valence.processor.fetch_materials_project_inputs", return_value=materials):
                results = processor.process_cation_system("Li", "O")

            self.assertEqual(results, [new_result])
            processor.service.fit_many.assert_called_once()
            pending_materials = processor.service.fit_many.call_args.args[0]
            self.assertEqual([material.material_id for material in pending_materials], ["mp-new"])

            dict_sijs = json.loads((pair_dir / "dict_sijs.json").read_text(encoding="utf-8"))
            dict_charges = json.loads((pair_dir / "dict_charges.json").read_text(encoding="utf-8"))
            self.assertEqual(sorted(dict_sijs), ["mp-new", "mp-solved"])
            self.assertEqual(sorted(dict_charges), ["mp-new", "mp-solved"])

            for algorithm in ("shgo", "diff"):
                self.assertTrue((pair_dir / "R0Bs" / algorithm / "mp-solved.txt").exists())
                self.assertTrue((pair_dir / "R0Bs" / algorithm / "mp-new.txt").exists())
                no_solution_text = (pair_dir / "no_solu" / f"{algorithm}.txt").read_text(encoding="utf-8")
                self.assertIn("mp-failed", no_solution_text)


if __name__ == "__main__":
    unittest.main()
