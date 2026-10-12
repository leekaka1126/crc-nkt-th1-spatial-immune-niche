"""Synthetic arithmetic checks for the saved-count E4 aggregation auditor."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from check_crc_R2k4_patient_estimand_independent import pair_effect  # noqa: E402


def row(passed: bool, edges: int, observed: int, shifted: int) -> dict[str, str]:
    return {"pair_passed_minimum": str(passed), "n_pair_retained": str(edges),
            "n_observed_SS_high_on_retained": str(observed),
            "n_shifted_SS_high_on_retained": str(shifted)}


class PairEffectTest(unittest.TestCase):
    def test_difference_of_differences_and_minimum_group_weight(self) -> None:
        errors: list[str] = []
        result = pair_effect({False: row(True, 40, 8, 4), True: row(True, 30, 12, 6)},
                             "synthetic", errors)
        self.assertEqual(errors, [])
        assert result is not None
        self.assertAlmostEqual(result[0], (12 / 30 - 8 / 40) - (6 / 30 - 4 / 40))
        self.assertEqual(result[1], 30)

    def test_failed_pair_is_missing_not_zero(self) -> None:
        errors: list[str] = []
        result = pair_effect({False: row(False, 0, 0, 0), True: row(False, 0, 0, 0)},
                             "synthetic", errors)
        self.assertIsNone(result)
        self.assertEqual(errors, [])

    def test_inconsistent_pair_gate_is_flagged(self) -> None:
        errors: list[str] = []
        pair_effect({False: row(False, 0, 0, 0), True: row(True, 30, 10, 5)},
                    "synthetic", errors)
        self.assertTrue(any("gate disagrees" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
