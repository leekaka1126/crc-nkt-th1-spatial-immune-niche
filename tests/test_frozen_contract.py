"""Small dependency-light checks of the locked panel, shifts, and patient unit."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def import_path(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


panel = import_path("panel_under_test", ROOT / "manuscript/CRC_NKT_v11_1/code_release/scripts/05_spatial_neighborhood_analysis.py")
shift = import_path("shift_under_test", ROOT / "scripts/gse294385_spatial_autocorr_audit.py")
e4 = import_path("e4_under_test", ROOT / "scripts/run_valdeolivas_pathology_matched_replication.py")


class FrozenContract(unittest.TestCase):
    def test_panels_remain_18_plus_13_unique(self):
        self.assertEqual(len(panel.IMMUNE_GENES), 18)
        self.assertEqual(len(panel.MATRIX_GENES), 13)
        self.assertEqual(len(set(panel.IMMUNE_GENES + panel.MATRIX_GENES)), 31)

    def test_offsets_are_fixed_and_balanced(self):
        offsets = shift.select_offsets()
        self.assertEqual(len(offsets), 96)
        self.assertEqual(len(set(offsets)), 96)
        self.assertEqual({radius: sum(row[2] == radius for row in offsets)
                          for radius in range(5, 11)}, {radius: 16 for radius in range(5, 11)})
        self.assertEqual(offsets, shift.select_offsets())

    def test_patient_equal_not_section_weighted(self):
        sections = pd.DataFrame({
            "dataset": ["synthetic"] * 3,
            "tissue_context": ["Primary_CRC"] * 3,
            "status": ["PASS"] * 3,
            "patient_id": ["P01", "P01", "P02"],
            "pathology_matched_shift_effect_fraction": [0.1, 0.3, 0.6],
        })
        summary, patients = e4.patient_summary(sections, "synthetic", "all")
        self.assertEqual(summary["n_eligible_sections"], 3)
        self.assertEqual(summary["n_eligible_patients"], 2)
        self.assertAlmostEqual(summary["estimate_fraction"], 0.4)
        self.assertEqual(len(patients), 2)


if __name__ == "__main__":
    unittest.main()
