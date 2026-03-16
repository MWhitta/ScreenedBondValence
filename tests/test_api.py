from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

from pymatgen.core import Lattice, Structure

from screened_bond_valence import build_material_from_cif, build_summary_payload
from screened_bond_valence.api import _prepare_structure_with_oxidation_states
from screened_bond_valence.models import (
    AlgorithmFitResult,
    BondValenceMaterial,
    MaterialFitResult,
    TheoreticalBondValenceResult,
)


def _make_result(material_id: str, formula: str) -> MaterialFitResult:
    material = BondValenceMaterial(
        material_id=material_id,
        cation="Li",
        anion="O",
        structure=object(),
        structure_graph=object(),
        formula_pretty=formula,
    )
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
        fits=(
            AlgorithmFitResult("shgo", 1.2, 0.45),
            AlgorithmFitResult("diff", 1.4, 0.55),
        ),
    )


def _make_lio_structure(*, with_oxidation_states: bool) -> Structure:
    species = ["Li+", "O2-"] if with_oxidation_states else ["Li", "O"]
    return Structure(
        Lattice.cubic(4.0),
        species,
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )


class BuildSummaryPayloadTests(unittest.TestCase):
    def test_supports_classifier_and_serializer(self) -> None:
        results = (
            _make_result("mp-oxide", "Li2O"),
            _make_result("mp-hydroxide", "LiOH"),
        )

        payload = build_summary_payload(
            results,
            classifier=lambda summary: "hydroxides" if "OH" in (summary.formula_pretty or "") else "oxides",
            serializer=lambda summary: summary.to_compact_dict(),
        )

        self.assertEqual(sorted(payload), ["hydroxides", "oxides"])
        oxide = payload["oxides"][0]
        self.assertEqual(oxide["mid"], "mp-oxide")
        self.assertEqual(oxide["formula"], "Li2O")
        self.assertAlmostEqual(oxide["R0"], 1.3)
        self.assertAlmostEqual(oxide["B"], 0.5)
        self.assertAlmostEqual(oxide["R0_std"], 0.1)
        self.assertAlmostEqual(oxide["B_std"], 0.05)
        self.assertEqual(oxide["n_algos"], 2)
        self.assertEqual(payload["hydroxides"][0]["mid"], "mp-hydroxide")


class BuildMaterialFromCifTests(unittest.TestCase):
    def test_builds_sample_material(self) -> None:
        cif_path = Path(__file__).resolve().parents[1] / "1011191_aspod.cif"
        material = build_material_from_cif(cif_path, cation="Li", anion="O")

        self.assertEqual(material.material_id, "1011191_aspod")
        self.assertEqual(material.cation, "Li")
        self.assertEqual(material.anion, "O")
        self.assertTrue(material.formula_pretty)
        self.assertGreater(len(material.structure), 0)
        self.assertGreater(len(material.structure_graph.graph.edges()), 0)
        self.assertIn(material.metadata["oxidation_state_source"], {"cif", "bv_analyzer", "composition_guess", "charge_map", "charge_map_override"})


class OxidationStatePreparationTests(unittest.TestCase):
    def test_prefers_bv_analyzer_for_missing_oxidation_states(self) -> None:
        undecorated = _make_lio_structure(with_oxidation_states=False)
        decorated = _make_lio_structure(with_oxidation_states=True)

        with patch(
            "screened_bond_valence.api.BVAnalyzer.get_oxi_state_decorated_structure",
            return_value=decorated,
        ) as mock_guess:
            prepared, charge_map, source = _prepare_structure_with_oxidation_states(undecorated)

        mock_guess.assert_called_once()
        self.assertEqual(source, "bv_analyzer")
        self.assertEqual(charge_map, {"Li": 1.0, "O": -2.0})
        self.assertEqual(prepared.composition.reduced_formula, undecorated.composition.reduced_formula)

    def test_falls_back_to_composition_guess(self) -> None:
        undecorated = Structure(
            Lattice.cubic(4.0),
            ["Li", "Li", "O"],
            [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]],
        )

        with patch(
            "screened_bond_valence.api.BVAnalyzer.get_oxi_state_decorated_structure",
            side_effect=ValueError("failed"),
        ):
            _, charge_map, source = _prepare_structure_with_oxidation_states(undecorated)

        self.assertEqual(source, "composition_guess")
        self.assertEqual(charge_map, {"Li": 1.0, "O": -2.0})


if __name__ == "__main__":
    unittest.main()
