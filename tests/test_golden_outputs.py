from __future__ import annotations

import json
import unittest
from pathlib import Path

from screened_bond_valence import ScreenedBondValenceService, build_material_from_cif, build_summary_payload


GOLDEN_DIR = Path(__file__).resolve().parent / "golden"
REPO_ROOT = Path(__file__).resolve().parents[1]


class GoldenOutputTests(unittest.TestCase):
    def _assert_matches_golden(self, cif_filename: str, golden_filename: str) -> None:
        golden = json.loads((GOLDEN_DIR / golden_filename).read_text(encoding="utf-8"))
        material = build_material_from_cif(REPO_ROOT / cif_filename, cation="Li", anion="O")
        result = ScreenedBondValenceService(algorithms=("shgo",)).fit_material(material)
        payload = build_summary_payload(
            [result],
            serializer=lambda summary: summary.to_compact_dict(),
        )

        self.assertEqual(material.material_id, golden["material_id"])
        self.assertEqual(material.formula_pretty, golden["formula_pretty"])
        self.assertEqual(material.metadata["oxidation_state_source"], golden["oxidation_state_source"])
        self.assertEqual(len(material.structure), golden["site_count"])
        self.assertEqual(len(material.structure_graph.graph.edges()), golden["edge_count"])

        self.assertIsNotNone(result.theoretical)
        assert result.theoretical is not None
        self.assertTrue(result.theoretical.has_solution)
        self.assertEqual(result.theoretical.charge_map, golden["charge_map"])
        self.assertEqual(result.theoretical.bond_valences, golden["bond_valences"])

        self.assertEqual(list(payload), ["results"])
        summary = payload["results"][0]
        expected_summary = golden["summary"]
        self.assertEqual(summary["mid"], expected_summary["mid"])
        self.assertEqual(summary["formula"], expected_summary["formula"])
        self.assertAlmostEqual(summary["R0"], expected_summary["R0"], places=8)
        self.assertAlmostEqual(summary["B"], expected_summary["B"], places=8)
        self.assertAlmostEqual(summary["R0_std"], expected_summary["R0_std"], places=8)
        self.assertAlmostEqual(summary["B_std"], expected_summary["B_std"], places=8)
        self.assertEqual(summary["n_algos"], expected_summary["n_algos"])

    def test_spodumene_sample_matches_golden(self) -> None:
        self._assert_matches_golden("1011191_aspod.cif", "1011191_aspod_shgo.json")

    def test_example_sample_matches_golden(self) -> None:
        self._assert_matches_golden("example.cif", "example_shgo.json")


if __name__ == "__main__":
    unittest.main()
