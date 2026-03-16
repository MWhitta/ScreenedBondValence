from __future__ import annotations

import unittest

from screened_bond_valence.solvers import TheoreticalBondValenceSolver


class TheoreticalBondValenceSolverTests(unittest.TestCase):
    def test_accepts_in_memory_charge_sources(self) -> None:
        solver = TheoreticalBondValenceSolver(
            element2charge={"Li": 1.0, "O": -2.0},
            species_by_material={"mp-test": ["Li+", "O2-"]},
        )

        self.assertEqual(solver.resolve_charge_map("mp-test"), {"Li": 1.0, "O": -2.0})
        self.assertEqual(
            solver.resolve_charge_map("mp-other", charge_map={"Li": 0.5, "O": -1.0}),
            {"Li": 0.5, "O": -1.0},
        )
        self.assertEqual(solver.resolve_charge_map("mp-fallback"), {"Li": 1.0, "O": -2.0})


if __name__ == "__main__":
    unittest.main()
