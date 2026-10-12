#!/usr/bin/env python
"""Patient-aware descriptive summaries of the frozen-effect edge-selection audit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
PAIRS = ("tumor->tumor", "tumor->stroma", "stroma->tumor", "stroma->stroma")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def edge2(table: pd.DataFrame, prefix: str) -> pd.Series:
    return table[[f"{prefix}_hop_{bin_name}" for bin_name in ("0", "1", "2")]].sum(axis=1)


def make_pair_contributions(pair_rows: pd.DataFrame, passing: pd.DataFrame) -> pd.DataFrame:
    pair_rows = pair_rows.merge(passing[["dataset", "sample", "shift_array_row", "shift_array_col"]],
                                on=["dataset", "sample", "shift_array_row", "shift_array_col"],
                                how="inner", validate="many_to_one")
    cells = pair_rows.pivot(index=["dataset", "patient_id", "sample", "shift_array_row", "shift_array_col",
                                   "pathology_pair"], columns="IS_high_source",
                            values=["n_pair_retained", "n_observed_SS_high_on_retained",
                                    "n_shifted_SS_high_on_retained"]).fillna(0)
    cells.columns = [f"{name}_{'high' if is_high else 'low'}" for name, is_high in cells.columns]
    cells = cells.reset_index()
    cells["pair_weight"] = cells[["n_pair_retained_low", "n_pair_retained_high"]].min(axis=1)
    weight_sum = cells.groupby(["dataset", "sample", "shift_array_row", "shift_array_col"]).pair_weight.transform("sum")
    cells["weight_share"] = cells.pair_weight / weight_sum
    cells["observed_pair_contrast"] = np.where(cells.pair_weight > 0,
        cells.n_observed_SS_high_on_retained_high / cells.n_pair_retained_high -
        cells.n_observed_SS_high_on_retained_low / cells.n_pair_retained_low, 0.0)
    cells["shifted_pair_contrast"] = np.where(cells.pair_weight > 0,
        cells.n_shifted_SS_high_on_retained_high / cells.n_pair_retained_high -
        cells.n_shifted_SS_high_on_retained_low / cells.n_pair_retained_low, 0.0)
    cells["pair_contribution"] = cells.weight_share * (
        cells.observed_pair_contrast - cells.shifted_pair_contrast)
    if not np.isclose(cells.groupby(["dataset", "sample", "shift_array_row", "shift_array_col"]).weight_share.sum(),
                      1.0).all():
        raise AssertionError("Pair weights do not sum to one")
    return cells


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--e4-dir", type=Path, required=True)
    args = parser.parse_args()
    audit_dir = args.audit_dir
    frozen_path = args.e4_dir / "Table_all_by_section.csv"
    shifts_path = audit_dir / "Table_shift_edge_selection.csv"
    pairs_path = audit_dir / "Table_shift_pair_and_edge_hops.csv"
    sections_path = audit_dir / "Table_section_baselines.csv"
    shifts = pd.read_csv(shifts_path)
    pairs = pd.read_csv(pairs_path)
    sections = pd.read_csv(sections_path)
    frozen = pd.read_csv(frozen_path)
    if len(shifts) != 1920 or len(pairs) != 15360 or len(sections) != 20:
        raise AssertionError("Audit record count mismatch")
    section_keys = ["dataset", "sample"]
    frozen = frozen[section_keys + ["status", "pathology_matched_shift_effect_fraction"]].rename(
        columns={"status": "frozen_section_status", "pathology_matched_shift_effect_fraction": "frozen_effect"})
    sections = sections.merge(frozen, on=section_keys, validate="one_to_one")
    if not sections.section_status.eq(sections.frozen_section_status).all():
        raise AssertionError("Frozen section status mismatch")
    accepted = shifts.loc[shifts.status.eq("PASS")].merge(
        sections[section_keys + ["section_status"]], on=section_keys, validate="many_to_one")
    accepted = accepted.loc[accepted.section_status.eq("PASS")].copy()
    shift_keys = ["dataset", "sample", "shift_array_row", "shift_array_col"]
    high_stages = pairs.loc[pairs.IS_high_source].groupby(shift_keys, as_index=False)[
        ["n_original_SS_low_edges", "n_mapped_in_tissue", "n_compatible", "n_pair_retained"]].sum()
    accepted = accepted.merge(high_stages, on=shift_keys, validate="one_to_one")
    if not accepted.n_original_SS_low_edges.eq(accepted.n_original_SS_low_IS_high_edges).all():
        raise AssertionError("High-edge original denominator mismatch")
    if not accepted.n_pair_retained.eq(accepted.n_pair_retained_IS_high_edges).all():
        raise AssertionError("High-edge retained denominator mismatch")
    accepted["mapped_high_fraction"] = accepted.n_mapped_in_tissue / accepted.n_original_SS_low_edges
    accepted["compatible_high_fraction"] = accepted.n_compatible / accepted.n_original_SS_low_edges
    accepted["retained_high_fraction"] = accepted.n_pair_retained / accepted.n_original_SS_low_edges
    accepted["off_tissue_high_loss_fraction"] = 1 - accepted.mapped_high_fraction
    accepted["pathology_mismatch_high_loss_fraction"] = (
        accepted.n_mapped_in_tissue - accepted.n_compatible) / accepted.n_original_SS_low_edges
    accepted["pair_minimum_high_loss_fraction"] = (
        accepted.n_compatible - accepted.n_pair_retained) / accepted.n_original_SS_low_edges
    accepted["effect_pp"] = 100 * accepted.observed_minus_shifted_IS_contrast
    accepted["observed_pp"] = 100 * accepted.observed_matched_edge_IS_contrast
    accepted["shifted_pp"] = 100 * accepted.shifted_matched_edge_IS_contrast
    accepted["observed_drift_pp"] = 100 * accepted.observed_matched_minus_full_edge_contrast
    accepted["original_high_edge2_fraction"] = edge2(accepted, "original_high_original_target") / (
        accepted.n_original_SS_low_IS_high_edges)
    accepted["retained_high_original_edge2_fraction"] = edge2(accepted, "retained_high_original_target") / (
        accepted.n_pair_retained_IS_high_edges)
    accepted["retained_high_mapped_edge2_fraction"] = edge2(accepted, "retained_high_mapped_target") / (
        accepted.n_pair_retained_IS_high_edges)
    accepted["original_low_edge2_fraction"] = edge2(accepted, "original_low_original_target") / (
        accepted.n_original_SS_low_source_edges - accepted.n_original_SS_low_IS_high_edges)
    accepted["retained_low_original_edge2_fraction"] = edge2(accepted, "retained_low_original_target") / (
        accepted.n_pair_retained_SS_low_source_edges - accepted.n_pair_retained_IS_high_edges)
    metrics = ["effect_pp", "observed_pp", "shifted_pp", "observed_drift_pp",
               "mapped_high_fraction", "compatible_high_fraction", "retained_high_fraction",
               "off_tissue_high_loss_fraction", "pathology_mismatch_high_loss_fraction",
               "pair_minimum_high_loss_fraction",
               "matched_IS_high_edge_fraction_after_pair_filter", "original_high_edge2_fraction",
               "retained_high_original_edge2_fraction", "retained_high_mapped_edge2_fraction",
               "original_low_edge2_fraction", "retained_low_original_edge2_fraction"]
    section_mean = accepted.groupby(["dataset", "patient_id", "sample"], as_index=False)[metrics].mean()
    section_count = accepted.groupby(section_keys).size().rename("n_accepted_shifts").reset_index()
    section_mean = sections.merge(section_mean, on=["dataset", "patient_id", "sample"], how="left",
                                  validate="one_to_one").merge(section_count, on=section_keys, how="left",
                                                                 suffixes=("_frozen", "_recount"), validate="one_to_one")
    good = section_mean.loc[section_mean.section_status.eq("PASS")]
    if not np.isclose(good.effect_pp / 100, good.frozen_effect, atol=1e-12).all():
        raise AssertionError("Accepted section mean differs from frozen effect")
    section_mean.to_csv(audit_dir / "Table_section_selection_summary.csv", index=False)
    patient_mean = good.groupby(["dataset", "patient_id"], as_index=False)[metrics].mean()
    section_number = good.groupby(["dataset", "patient_id"]).size().rename("n_passing_sections").reset_index()
    patient_mean = patient_mean.merge(section_number, on=["dataset", "patient_id"], validate="one_to_one")
    patient_mean.to_csv(audit_dir / "Table_patient_selection_summary.csv", index=False)
    cohort_mean = patient_mean.groupby("dataset", as_index=False)[metrics].mean()
    cohort_mean["n_passing_patients"] = patient_mean.groupby("dataset").size().to_numpy()
    cohort_mean["n_passing_sections"] = good.groupby("dataset").size().to_numpy()
    cohort_mean.to_csv(audit_dir / "Table_patient_equal_cohort_selection.csv", index=False)
    cells = make_pair_contributions(pairs, accepted)
    section_pairs = cells.groupby(["dataset", "patient_id", "sample", "pathology_pair"], as_index=False)[
        ["weight_share", "pair_contribution", "n_pair_retained_high", "n_pair_retained_low"]].mean()
    section_pairs.to_csv(audit_dir / "Table_pair_contribution_by_section.csv", index=False)
    patient_pairs = section_pairs.groupby(["dataset", "patient_id", "pathology_pair"], as_index=False)[
        ["weight_share", "pair_contribution"]].mean()
    patient_pairs.to_csv(audit_dir / "Table_pair_contribution_by_patient.csv", index=False)
    check_contrib = patient_pairs.groupby("dataset").pair_contribution.mean() * len(PAIRS)
    check_effect = patient_mean.groupby("dataset").effect_pp.mean() / 100
    if not np.allclose(check_contrib.sort_index(), check_effect.sort_index(), atol=1e-12):
        raise AssertionError("Pair contributions do not sum to patient-equal effect")
    half_rows = []
    for (dataset, patient_id, sample), group in accepted.groupby(["dataset", "patient_id", "sample"]):
        group = group.sort_values(["matched_IS_high_edge_fraction_after_pair_filter",
                                   "shift_array_row", "shift_array_col"])
        n = len(group)
        lower, upper = group.iloc[:n // 2], group.iloc[n // 2:]
        rho_coverage = spearmanr(group.matched_IS_high_edge_fraction_after_pair_filter,
                                 group.effect_pp).statistic if n >= 3 else np.nan
        rho_edge = spearmanr(group.retained_high_original_edge2_fraction,
                             group.effect_pp).statistic if n >= 3 else np.nan
        half_rows.append({"dataset": dataset, "patient_id": patient_id, "sample": sample,
                          "n_accepted_shifts": n, "n_lower_coverage": len(lower),
                          "n_upper_coverage": len(upper),
                          "lower_coverage_mean_effect_pp": lower.effect_pp.mean(),
                          "upper_coverage_mean_effect_pp": upper.effect_pp.mean(),
                          "upper_minus_lower_effect_pp": upper.effect_pp.mean() - lower.effect_pp.mean(),
                          "spearman_coverage_effect": rho_coverage,
                          "spearman_edge2_effect": rho_edge})
    halves = pd.DataFrame(half_rows)
    halves.to_csv(audit_dir / "Table_section_coverage_halves.csv", index=False)
    halves.groupby(["dataset", "patient_id"], as_index=False)[
        ["lower_coverage_mean_effect_pp", "upper_coverage_mean_effect_pp",
         "upper_minus_lower_effect_pp", "spearman_coverage_effect",
         "spearman_edge2_effect"]].mean().to_csv(audit_dir / "Table_patient_coverage_halves.csv", index=False)
    manifest = {"analysis_role": "post_hoc_descriptive", "source_shift_sha256": sha256(shifts_path),
                "source_pair_sha256": sha256(pairs_path), "source_section_sha256": sha256(sections_path),
                "frozen_section_sha256": sha256(frozen_path), "script_sha256": sha256(Path(__file__)),
                "n_passing_shifts_in_passing_sections": len(accepted),
                "unit_of_cohort_summary": "patient, after section mean",
                "coverage_halves": "within-section lower and upper rank halves of accepted shifts; dependent post hoc diagnostic",
                "no_P_values_or_exact_permutation_claim": True}
    (audit_dir / "summary_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(cohort_mean.to_string(index=False))
    print("Cross-pair retained qualified edges in external passing sections:",
          int(cells.loc[cells.dataset.eq("Valdeolivas_2024") &
                        cells.pathology_pair.isin(["tumor->stroma", "stroma->tumor"]),
                        "pair_weight"].sum()))


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
