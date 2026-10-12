#!/usr/bin/env python
"""Audit pair estimability at section, patient and cohort levels without changing E4.

All identifiers are dataset-local public sample/patient codes. Offset counts are
dependent technical support, never an inferential sample size.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from summarize_v11_8_full_edge_selection import make_pair_contributions, PAIRS

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, required=True,
                        help="New edge-support directory produced by this public-input run")
    args = parser.parse_args()
    audit = args.audit_dir
    shift_path = audit / "Table_shift_edge_selection.csv"
    pair_path = audit / "Table_shift_pair_and_edge_hops.csv"
    section_path = audit / "Table_section_selection_summary.csv"
    patient_contrib_path = audit / "Table_pair_contribution_by_patient.csv"
    out_path = audit / "Table_S14i_multilevel_pair_estimability.csv"
    shifts = pd.read_csv(shift_path)
    pairs = pd.read_csv(pair_path)
    sections = pd.read_csv(section_path)
    recomputed_patient_pairs = pd.read_csv(patient_contrib_path)
    assert len(shifts) == 1920 and len(pairs) == 15360 and len(sections) == 20

    section_keys = ["dataset", "patient_id", "sample"]
    offset_keys = section_keys + ["shift_array_row", "shift_array_col"]
    passed_sections = sections.loc[sections.section_status.eq("PASS"), ["dataset", "sample"]]
    analysis_offsets = shifts.loc[shifts.status.eq("PASS")].merge(
        passed_sections, on=["dataset", "sample"], validate="many_to_one")
    assert len(analysis_offsets) == 1395  # 672 source + 723 independent
    cells = make_pair_contributions(pairs, analysis_offsets)
    cells["pair_specific_effect_pp"] = 100 * (
        cells.observed_pair_contrast - cells.shifted_pair_contrast)
    cells["mixture_contribution_pp"] = 100 * cells.pair_contribution

    # S14b has two IS-group rows per pair/offset. Count a qualified pair once.
    technical = pairs.loc[pairs.status.eq("PASS") & pairs.pair_passed_minimum,
                          offset_keys + ["pathology_pair"]].drop_duplicates()
    assert technical.duplicated(offset_keys + ["pathology_pair"]).sum() == 0
    records: list[dict] = []
    for section in sections.itertuples(index=False):
        key = (section.dataset, section.patient_id, section.sample)
        section_shifts = shifts.loc[(shifts.dataset == key[0]) &
                                    (shifts.patient_id == key[1]) & (shifts["sample"] == key[2])]
        n_accepted = int(section_shifts.status.eq("PASS").sum())
        assert len(section_shifts) == 96
        for pair in PAIRS:
            sub = cells.loc[(cells.dataset == key[0]) & (cells.patient_id == key[1]) &
                            (cells["sample"] == key[2]) & (cells.pathology_pair == pair)]
            qualified = sub.loc[sub.pair_weight.gt(0)]
            technical_n = len(technical.loc[(technical.dataset == key[0]) &
                                              (technical.patient_id == key[1]) &
                                              (technical["sample"] == key[2]) &
                                              (technical.pathology_pair == pair)])
            is_pass = section.section_status == "PASS"
            assert (len(sub) == n_accepted) if is_pass else sub.empty
            records.append({
                "scope": "section", "dataset": key[0], "patient_id": key[1],
                "sample": key[2], "pathology_pair": pair,
                "section_status": section.section_status,
                "attempted_offset_n": len(section_shifts),
                "offset_pass_n_before_section_gate": n_accepted,
                "analysis_offset_n": n_accepted if is_pass else 0,
                "technical_qualified_offset_n_before_section_gate": technical_n,
                "qualified_analysis_offset_n": len(qualified),
                "mean_retained_low_edges_given_qualified": qualified.n_pair_retained_low.mean(),
                "mean_retained_high_edges_given_qualified": qualified.n_pair_retained_high.mean(),
                "mean_effective_pair_weight_given_qualified": qualified.pair_weight.mean(),
                "mean_weight_share_over_analysis_offsets": sub.weight_share.mean(),
                "mean_mixture_contribution_pp": sub.mixture_contribution_pp.mean(),
                "pair_specific_effect_pp_given_qualified":
                    qualified.pair_specific_effect_pp.mean(),
                "pair_effect_status": ("SECTION_FAILED" if not is_pass else
                                       "NOT_ESTIMABLE" if qualified.empty else "ESTIMABLE"),
                "nonestimability_reason": ("Section failed the fixed >=24 accepted-offset gate"
                    if not is_pass else "No accepted offset met the >=30-edge minimum in both IS groups"
                    if qualified.empty else ""),
            })
    section_table = pd.DataFrame(records)

    # Preserve the original hierarchy: offsets -> passing sections -> patients
    # -> cohort. Conditional pair-specific effects are audit descriptors only.
    for scope, source, keys in (
        ("patient", section_table.loc[section_table.section_status.eq("PASS")],
         ["dataset", "patient_id", "pathology_pair"]),
        ("cohort", None, ["dataset", "pathology_pair"]),
    ):
        if source is None:
            source = pd.DataFrame(records).loc[lambda x: x.scope.eq("patient")]
        for group_key, group in source.groupby(keys, dropna=False):
            ident = dict(zip(keys, group_key))
            qualified = group.loc[group.qualified_analysis_offset_n.gt(0)]
            records.append({
                "scope": scope, "dataset": ident["dataset"],
                "patient_id": ident.get("patient_id", ""), "sample": "",
                "pathology_pair": ident["pathology_pair"], "section_status": "PASS",
                "attempted_offset_n": int(group.attempted_offset_n.sum()),
                "offset_pass_n_before_section_gate": int(group.offset_pass_n_before_section_gate.sum()),
                "analysis_offset_n": int(group.analysis_offset_n.sum()),
                "technical_qualified_offset_n_before_section_gate":
                    int(group.technical_qualified_offset_n_before_section_gate.sum()),
                "qualified_analysis_offset_n": int(group.qualified_analysis_offset_n.sum()),
                "mean_retained_low_edges_given_qualified":
                    np.average(qualified.mean_retained_low_edges_given_qualified,
                               weights=qualified.qualified_analysis_offset_n) if len(qualified) else np.nan,
                "mean_retained_high_edges_given_qualified":
                    np.average(qualified.mean_retained_high_edges_given_qualified,
                               weights=qualified.qualified_analysis_offset_n) if len(qualified) else np.nan,
                "mean_effective_pair_weight_given_qualified":
                    np.average(qualified.mean_effective_pair_weight_given_qualified,
                               weights=qualified.qualified_analysis_offset_n) if len(qualified) else np.nan,
                "mean_weight_share_over_analysis_offsets": group.mean_weight_share_over_analysis_offsets.mean(),
                "mean_mixture_contribution_pp": group.mean_mixture_contribution_pp.mean(),
                "pair_specific_effect_pp_given_qualified":
                    qualified.pair_specific_effect_pp_given_qualified.mean() if len(qualified) else np.nan,
                "pair_effect_status": "ESTIMABLE" if len(qualified) else "NOT_ESTIMABLE",
                "nonestimability_reason": "" if len(qualified) else
                    "No qualifying pair among contributing offsets",
            })
    out = pd.DataFrame(records)
    # A target contrast can be mathematically specified yet unavailable in the
    # observed data. Do not conflate estimand definition with pair eligibility.
    out["estimand_defined"] = True
    out["pair_eligible"] = out.pair_effect_status.eq("ESTIMABLE")
    out["estimate"] = out.pair_specific_effect_pp_given_qualified
    out["mixture_weight"] = out.mean_weight_share_over_analysis_offsets
    out["contribution"] = out.mean_mixture_contribution_pp
    out["not_estimable_reason"] = out.nonestimability_reason
    front = ["scope", "dataset", "patient_id", "sample", "pathology_pair",
             "section_status", "estimand_defined", "pair_eligible", "estimate",
             "mixture_weight", "contribution", "not_estimable_reason"]
    out = out[front + [col for col in out if col not in front]]
    assert len(out) == 80 + 44 + 8
    assert out.loc[(out.dataset == "Valdeolivas_2024") &
                   out.pathology_pair.isin(["tumor->stroma", "stroma->tumor"]) &
                   out.scope.eq("cohort"), "pair_effect_status"].eq("NOT_ESTIMABLE").all()
    external_cross = out.loc[out.scope.eq("cohort") & out.dataset.eq("Valdeolivas_2024") &
                             out.pathology_pair.isin(["tumor->stroma", "stroma->tumor"])]
    assert external_cross.estimand_defined.all() and not external_cross.pair_eligible.any()
    assert external_cross.estimate.isna().all()
    assert external_cross.mixture_weight.eq(0).all() and external_cross.contribution.eq(0).all()
    patient = out.loc[out.scope.eq("patient")].merge(
        recomputed_patient_pairs, on=["dataset", "patient_id", "pathology_pair"],
        validate="one_to_one")
    assert np.allclose(patient.mean_mixture_contribution_pp,
                       100 * patient.pair_contribution, atol=1e-10)
    assert np.allclose(patient.mean_weight_share_over_analysis_offsets,
                       patient.weight_share, atol=1e-10)
    cohort_effect = out.loc[out.scope.eq("cohort")].groupby("dataset").mean_mixture_contribution_pp.sum()
    recomputed_effect = sections.loc[sections.section_status.eq("PASS")].groupby(
        ["dataset", "patient_id"]).effect_pp.mean().groupby("dataset").mean()
    assert np.allclose(cohort_effect.sort_index(), recomputed_effect.sort_index(), atol=1e-10)
    out.to_csv(out_path, index=False, float_format="%.12g")
    manifest = {"role": "post_hoc_pair_support_audit_not_new_primary_estimate",
                "output_sha256": sha256(out_path), "script_sha256": sha256(Path(__file__)),
                "source_sha256": {
                    p.name: sha256(p) for p in
                    (shift_path, pair_path, section_path, patient_contrib_path)},
                "rows": len(out), "section_rows": 80, "patient_rows": 44, "cohort_rows": 8,
                "analysis_offsets": len(analysis_offsets),
                "nonestimability": "NA pair-specific effect; zero mixture contribution only means pair excluded",
                "estimand_definition": "The theoretical pair contrast is specified for all four pairs; eligibility, not mathematical definition, fails when support is absent",
                "aggregation": "offset mean within passing section, equal section within patient, equal patient within cohort"}
    (audit / "Table_S14i_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
