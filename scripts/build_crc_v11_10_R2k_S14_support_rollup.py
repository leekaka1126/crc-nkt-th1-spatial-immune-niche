#!/usr/bin/env python
"""Analysis: fixed post hoc S14b pre/post-gate support rollup.

Date: 2026-10-11
Random seed: 42 (deterministic counts; no random operation is used)
Python and package versions: recorded in the output manifest at execution.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


np.random.seed(42)
OFFSET_KEYS = ["shift_array_row", "shift_array_col"]
PAIR_ORDER = ["tumor->tumor", "tumor->stroma", "stroma->tumor", "stroma->stroma"]
GATE_COLS = ["n_mapped_in_tissue", "n_off_tissue_or_lattice",
             "n_pathology_mismatch", "n_compatible"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_input(b: pd.DataFrame) -> None:
    required = ["dataset", "sample", "patient_id", *OFFSET_KEYS, "status", "pathology_pair",
                "IS_high_source", "n_original_SS_low_edges", *GATE_COLS,
                "pair_passed_minimum", "n_pair_retained"]
    missing = [column for column in required if column not in b]
    if missing or b[required].isna().any().any():
        raise AssertionError(f"Missing required fields or values: {missing}")
    row_key = ["dataset", "sample", *OFFSET_KEYS, "pathology_pair", "IS_high_source"]
    if len(b) != 15360 or b.duplicated(row_key).any():
        raise AssertionError("Expected 15,360 unique section-offset-pair-group rows")
    if set(b.pathology_pair) != set(PAIR_ORDER):
        raise AssertionError("Unexpected directed pathology-pair set")
    section_offset = b.groupby(["dataset", "sample"], sort=False)[OFFSET_KEYS[0]].size()
    if len(section_offset) != 20 or not section_offset.eq(96 * 4 * 2).all():
        raise AssertionError("Each of 20 sections must have 96 offsets × four pairs × two groups")
    offsets = b.groupby(["dataset", "sample", *OFFSET_KEYS], sort=False)
    if len(offsets) != 1920 or not offsets.size().eq(8).all() or not offsets.status.nunique().eq(1).all():
        raise AssertionError("Each offset must have eight pair/group cells and one frozen status")
    if (b[["n_original_SS_low_edges", *GATE_COLS, "n_pair_retained"]] < 0).any().any():
        raise AssertionError("Negative edge count")
    if not (b.n_original_SS_low_edges == b.n_mapped_in_tissue + b.n_off_tissue_or_lattice).all():
        raise AssertionError("Original-to-mapped edge flow does not partition")
    if not (b.n_mapped_in_tissue == b.n_compatible + b.n_pathology_mismatch).all():
        raise AssertionError("Mapped-to-compatible edge flow does not partition")
    if not (b.n_pair_retained <= b.n_compatible).all():
        raise AssertionError("Retained edges exceed compatible edges")
    if (b.loc[~b.pair_passed_minimum, "n_pair_retained"] != 0).any():
        raise AssertionError("Pair-failed row retains edges")


def section_rollup(b: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, sample, patient_id, pair), frame in b.groupby(
            ["dataset", "sample", "patient_id", "pathology_pair"], sort=True):
        low = frame.loc[~frame.IS_high_source].set_index(OFFSET_KEYS).sort_index()
        high = frame.loc[frame.IS_high_source].set_index(OFFSET_KEYS).sort_index()
        if len(low) != 96 or len(high) != 96 or not low.index.equals(high.index):
            raise AssertionError(f"Unmatched offsets or source-IS groups: {dataset}/{sample}/{pair}")
        if not low.status.equals(high.status) or not low.pair_passed_minimum.equals(high.pair_passed_minimum):
            raise AssertionError(f"Inconsistent status or pair gate between groups: {dataset}/{sample}/{pair}")
        passed = low.status.eq("PASS")
        pair_passed = low.pair_passed_minimum & high.pair_passed_minimum
        qualified_before_section = passed & pair_passed
        offset_pass_n = int(passed.sum())
        section_pass = offset_pass_n >= 24
        qualified_analysis = qualified_before_section if section_pass else pd.Series(False, index=low.index)
        row = {
            "dataset": dataset, "sample": sample, "patient_id": patient_id,
            "pathology_pair": pair, "section_status": "PASS" if section_pass else "FAIL",
            "attempted_offset_n": 96,
            "offset_pass_before_section_gate_n": offset_pass_n,
            "pair_qualified_before_section_gate_n": int(qualified_before_section.sum()),
            "analysis_offset_n": offset_pass_n if section_pass else 0,
            "pair_qualified_analysis_offset_n": int(qualified_analysis.sum()),
            "pair_support_status": ("SECTION_FAILED" if not section_pass else
                                    "QUALIFIED" if qualified_analysis.any() else "NOT_ESTIMABLE"),
        }
        for label, group in (("low", low), ("high", high)):
            original = group.n_original_SS_low_edges.unique()
            if len(original) != 1:
                raise AssertionError(f"Original edge count varies across offsets: {dataset}/{sample}/{pair}/{label}")
            row[f"original_distinct_edges_{label}"] = int(original[0])
            for column in GATE_COLS:
                row[f"median_{column}_per_attempted_offset_{label}"] = float(group[column].median())
            if section_pass:
                row[f"retained_edge_offset_incidence_in_analysis_{label}"] = int(
                    group.loc[qualified_analysis, "n_pair_retained"].sum()
                )
            else:
                row[f"retained_edge_offset_incidence_in_analysis_{label}"] = np.nan
            row[f"median_retained_edges_per_qualified_analysis_offset_{label}"] = (
                float(group.loc[qualified_analysis, "n_pair_retained"].median())
                if qualified_analysis.any() else np.nan
            )
        rows.append(row)
    section = pd.DataFrame(rows)
    section["pathology_pair"] = pd.Categorical(section.pathology_pair, PAIR_ORDER, ordered=True)
    section = section.sort_values(["dataset", "patient_id", "sample", "pathology_pair"])
    section["pathology_pair"] = section.pathology_pair.astype(str)
    return section.reset_index(drop=True)


def patient_rollup(section: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, patient_id, pair), frame in section.groupby(
            ["dataset", "patient_id", "pathology_pair"], sort=True):
        passing = frame.loc[frame.section_status.eq("PASS")]
        qualified = passing.loc[passing.pair_qualified_analysis_offset_n.gt(0)]
        row = {
            "dataset": dataset, "patient_id": patient_id, "pathology_pair": pair,
            "candidate_section_n": len(frame), "passing_section_n": len(passing),
            "sections_with_qualified_pair_n": len(qualified),
            "candidate_attempted_offset_n": int(frame.attempted_offset_n.sum()),
            "offset_pass_before_section_gate_all_candidates_n": int(
                frame.offset_pass_before_section_gate_n.sum()),
            "pair_qualified_before_section_gate_all_candidates_n": int(
                frame.pair_qualified_before_section_gate_n.sum()),
            "analysis_offset_n": int(frame.analysis_offset_n.sum()),
            "pair_qualified_analysis_offset_n": int(frame.pair_qualified_analysis_offset_n.sum()),
            "pair_support_status": ("NO_PASSING_SECTION" if passing.empty else
                                    "QUALIFIED" if not qualified.empty else "NOT_ESTIMABLE"),
        }
        for label in ("low", "high"):
            row[f"original_distinct_edges_across_sections_{label}"] = int(
                frame[f"original_distinct_edges_{label}"].sum()
            )
            row[f"retained_edge_offset_incidence_in_analysis_{label}"] = (
                int(passing[f"retained_edge_offset_incidence_in_analysis_{label}"].sum())
                if not passing.empty else np.nan
            )
            row[f"median_of_section_medians_retained_edges_given_qualified_{label}"] = (
                float(qualified[f"median_retained_edges_per_qualified_analysis_offset_{label}"].median())
                if not qualified.empty else np.nan
            )
        rows.append(row)
    patient = pd.DataFrame(rows)
    patient["pathology_pair"] = pd.Categorical(patient.pathology_pair, PAIR_ORDER, ordered=True)
    patient = patient.sort_values(["dataset", "patient_id", "pathology_pair"])
    patient["pathology_pair"] = patient.pathology_pair.astype(str)
    return patient.reset_index(drop=True)


def reconcile(section: pd.DataFrame, patient: pd.DataFrame) -> dict:
    """Check denominator/status arithmetic without reading any archived ledger."""
    if len(section) != 80 or len(patient) != 44:
        raise AssertionError("Unexpected section/patient/pair denominator")
    totals = section.groupby(["dataset", "patient_id", "pathology_pair"], as_index=False)[
        ["attempted_offset_n", "analysis_offset_n", "pair_qualified_analysis_offset_n"]].sum()
    merged = patient.merge(totals, on=["dataset", "patient_id", "pathology_pair"],
                           validate="one_to_one", suffixes=("_patient", "_section"))
    if len(merged) != len(patient):
        raise AssertionError("Patient rollup lost section keys")
    for field in ("analysis_offset_n", "pair_qualified_analysis_offset_n"):
        if not merged[f"{field}_patient"].eq(merged[f"{field}_section"]).all():
            raise AssertionError(f"Patient/section total mismatch: {field}")
    if not patient.candidate_attempted_offset_n.eq(
            merged.attempted_offset_n).all():
        raise AssertionError("Candidate attempt total mismatch")
    cohort = patient.groupby(["dataset", "pathology_pair"], sort=True).agg(
        analysis_offset_n=("analysis_offset_n", "sum"),
        pair_qualified_analysis_offset_n=("pair_qualified_analysis_offset_n", "sum"),
        patient_n=("patient_id", "nunique"),
    ).reset_index()
    for dataset, denominator, patient_n in (("GSE294385", 672, 5), ("Valdeolivas_2024", 723, 6)):
        rows = cohort.loc[cohort.dataset.eq(dataset)]
        if len(rows) != 4 or not rows.analysis_offset_n.eq(denominator).all() or not rows.patient_n.eq(patient_n).all():
            raise AssertionError(f"Cohort denominator mismatch: {dataset}")
    cross = cohort.loc[cohort.dataset.eq("Valdeolivas_2024") &
                       cohort.pathology_pair.isin(["tumor->stroma", "stroma->tumor"])]
    if len(cross) != 2 or not cross.pair_qualified_analysis_offset_n.eq(0).all():
        raise AssertionError("Independent cross-compartment non-estimability changed")
    passing_original = section.loc[section.section_status.eq("PASS")].groupby(
        ["dataset", "pathology_pair"], sort=True
    )[["original_distinct_edges_low", "original_distinct_edges_high"]].sum().reset_index()
    return {"section_rows": len(section), "patient_rows": len(patient),
            "cohort_pair_checks": cohort.to_dict(orient="records"),
            "original_edge_counts_across_passing_sections": passing_original.to_dict(orient="records")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, required=True)
    args = parser.parse_args()
    source = args.audit_dir / "Table_shift_pair_and_edge_hops.csv"
    section_out = args.audit_dir / "Table_S14p_pre_post_gate_support_by_section.csv"
    patient_out = args.audit_dir / "Table_S14q_pre_post_gate_support_by_patient.csv"
    manifest_out = args.audit_dir / "Table_S14p-q_support_rollup_manifest.json"
    b = pd.read_csv(source)
    check_input(b)
    section = section_rollup(b)
    patient = patient_rollup(section)
    summary = reconcile(section, patient)
    section.to_csv(section_out, index=False, float_format="%.3f")
    patient.to_csv(patient_out, index=False, float_format="%.3f")
    manifest = {
        "status": "PASS_SUPPORT_ONLY_NO_EFFECT_REESTIMATION",
        "plan": "post hoc support-only rollup of the current raw-input E4 edge ledger",
        "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
        "inputs": {source.name: sha256(source)},
        "outputs": {path.name: sha256(path) for path in (section_out, patient_out)},
        "summary": summary,
        "interpretation": "Original edges counted once per section; retained edge-offset incidence repeats edges across technical offsets; no pair effect estimated."
    }
    manifest_out.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"PASS: {len(section)} section × pair rows; {len(patient)} patient × pair rows")
    print(cohort_markdown(patient, section))


def cohort_markdown(patient: pd.DataFrame, section: pd.DataFrame) -> str:
    cohort = patient.groupby(["dataset", "pathology_pair"], sort=True).agg(
        patients=("patient_id", "nunique"),
        analysis_offsets=("analysis_offset_n", "sum"),
        qualified_offsets=("pair_qualified_analysis_offset_n", "sum"),
    ).reset_index()
    pre = section.loc[section.section_status.eq("PASS")].groupby(
        ["dataset", "pathology_pair"], sort=True
    )[["original_distinct_edges_low", "original_distinct_edges_high"]].sum().reset_index()
    cohort = cohort.merge(pre, on=["dataset", "pathology_pair"], validate="one_to_one")
    lines = ["| dataset | pathology pair | patients | original IS-low / IS-high edges in passing sections | analysis offsets | qualified offsets |",
             "|---|---|---:|---:|---:|---:|"]
    for row in cohort.itertuples(index=False):
        lines.append(f"| {row.dataset} | {row.pathology_pair} | {row.patients} | "
                     f"{row.original_distinct_edges_low} / {row.original_distinct_edges_high} | "
                     f"{row.analysis_offsets} | {row.qualified_offsets} |")
    return "\n".join(lines)


if __name__ == "__main__":
    main()
