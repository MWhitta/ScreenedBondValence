from __future__ import annotations

import unittest

from screened_bond_valence.models import (
    AlgorithmFitResult,
    BondValenceMaterial,
    MaterialFitResult,
    TheoreticalBondValenceResult,
)


class MaterialFitSummaryTests(unittest.TestCase):
    def test_aggregate_and_compact_dict(self) -> None:
        material = BondValenceMaterial(
            material_id="mp-test",
            cation="Li",
            anion="O",
            structure=object(),
            structure_graph=object(),
            formula_pretty="Li2O",
            metadata={"source": "unit-test"},
        )
        theoretical = TheoreticalBondValenceResult(
            material_id="mp-test",
            bond_valences={"Li1O1": 1.0},
            bond_types=("Li1O1",),
            bond_lengths={"Li1O1": 2.0},
            charge_map={"Li": 1.0, "O": -2.0},
        )
        result = MaterialFitResult(
            material=material,
            theoretical=theoretical,
            fits=(
                AlgorithmFitResult("shgo", 1.1, 0.4),
                AlgorithmFitResult("diff", 1.5, 0.6),
                AlgorithmFitResult("direct", 1.3, 0.5),
            ),
        )

        summary = result.aggregate()

        self.assertIsNotNone(summary)
        assert summary is not None
        self.assertEqual(summary.material_id, "mp-test")
        self.assertAlmostEqual(summary.r0, 1.3)
        self.assertAlmostEqual(summary.b, 0.5)
        self.assertEqual(summary.n_algos, 3)
        self.assertEqual(summary.algorithms, ("shgo", "diff", "direct"))
        self.assertEqual(summary.metadata, {"source": "unit-test"})
        self.assertEqual(
            summary.to_compact_dict(),
            {
                "mid": "mp-test",
                "formula": "Li2O",
                "R0": 1.3,
                "B": 0.5,
                "R0_std": summary.r0_std,
                "B_std": summary.b_std,
                "n_algos": 3,
            },
        )


if __name__ == "__main__":
    unittest.main()
