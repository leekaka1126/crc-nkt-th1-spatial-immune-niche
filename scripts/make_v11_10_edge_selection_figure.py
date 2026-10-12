#!/usr/bin/env python
"""Create the v11.10 main edge-selection and pathology-pair support figure."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OKABE = {"source": "#0072B2", "external": "#D55E00", "gray": "#666666", "black": "#000000"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    args = parser.parse_args()
    audit_dir, figure_dir = args.audit_dir, args.figure_dir
    cohort = pd.read_csv(audit_dir / "Table_patient_equal_cohort_selection.csv")
    patients = pd.read_csv(audit_dir / "Table_patient_selection_summary.csv")
    halves = pd.read_csv(audit_dir / "Table_patient_coverage_halves.csv")
    shift_rows = pd.read_csv(audit_dir / "Table_shift_edge_selection.csv")
    pair_rows = pd.read_csv(audit_dir / "Table_shift_pair_and_edge_hops.csv")
    sections = pd.read_csv(audit_dir / "Table_section_selection_summary.csv")
    cohort = cohort.set_index("dataset").loc[["GSE294385", "Valdeolivas_2024"]]
    labels = {"GSE294385": "GSE294385 primary", "Valdeolivas_2024": "Valdeolivas"}
    colors = {"GSE294385": OKABE["source"], "Valdeolivas_2024": OKABE["external"]}
    figure_dir.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42,
                         "ps.fonttype": 42})
    fig, grid = plt.subplots(2, 2, figsize=(7.0, 5.5), gridspec_kw={"wspace": 0.42, "hspace": 0.60})
    axes = grid.ravel()
    stages = ["mapped_high_fraction", "compatible_high_fraction", "retained_high_fraction"]
    stage_labels = ["In-tissue\nspot", "Pathology\nmatch", "Pair\nminimum"]
    x = np.arange(3)
    width = 0.34
    for offset, dataset in zip((-width / 2, width / 2), ["GSE294385", "Valdeolivas_2024"]):
        values = 100 * cohort.loc[dataset, stages].to_numpy(float)
        axes[0].bar(x + offset, values, width, label=labels[dataset], color=colors[dataset],
                    edgecolor="white", linewidth=0.35)
    axes[0].set_xticks(x, stage_labels)
    axes[0].set_ylim(0, 100)
    axes[0].set_ylabel("Remaining IS-high edges (%)")
    axes[0].set_title("Edge flow", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=6.5, loc="upper right")
    axes[0].text(0.02, -0.28, "Patient-equal means; original SS-low-source IS-high denominator",
                 transform=axes[0].transAxes, fontsize=6.3, color=OKABE["gray"])

    for dataset in ["GSE294385", "Valdeolivas_2024"]:
        sub = patients.loc[patients.dataset.eq(dataset)]
        axes[1].scatter(100 * sub.retained_high_original_edge2_fraction, sub.effect_pp,
                        s=28, color=colors[dataset], edgecolor="white", linewidth=0.4,
                        label=labels[dataset], alpha=0.95)
        for row in sub.itertuples(index=False):
            if dataset == "Valdeolivas_2024" and row.patient_id == "A120838":
                axes[1].annotate("A120838", (100 * row.retained_high_original_edge2_fraction, row.effect_pp),
                                 xytext=(4, 4), textcoords="offset points", fontsize=6.4,
                                 color=colors[dataset])
    axes[1].axhline(0, color=OKABE["black"], linewidth=0.7)
    axes[1].set_xlabel("Retained targets within ≤2 edge hops (%)")
    axes[1].set_ylabel("Patient effect (percentage points)")
    axes[1].set_title("Edge location", loc="left", fontweight="bold")
    axes[1].set_xlim(0, max(60, float((100 * patients.retained_high_original_edge2_fraction).max() + 4)))
    axes[1].legend(frameon=False, fontsize=6.5, loc="lower right")

    # Coverage-half comparison is a dependent-shift diagnostic, not an uncertainty bar.
    halves = halves.copy()
    positions = {"GSE294385": 0, "Valdeolivas_2024": 1}
    for dataset in ["GSE294385", "Valdeolivas_2024"]:
        sub = halves.loc[halves.dataset.eq(dataset)].copy()
        xpos = positions[dataset]
        rng = np.random.default_rng(42 + xpos)
        jitter = rng.uniform(-0.07, 0.07, size=len(sub))
        axes[2].plot(np.full(len(sub), xpos - 0.13) + jitter, sub.lower_coverage_mean_effect_pp,
                     "o", ms=4.1, color=OKABE["gray"], alpha=0.75, label="Lower coverage" if xpos == 0 else None)
        axes[2].plot(np.full(len(sub), xpos + 0.13) + jitter, sub.upper_coverage_mean_effect_pp,
                     "o", ms=4.1, color=colors[dataset], alpha=0.85, label="Higher coverage" if xpos == 0 else None)
        for low, high, j in zip(sub.lower_coverage_mean_effect_pp, sub.upper_coverage_mean_effect_pp, jitter):
            axes[2].plot([xpos - 0.13 + j, xpos + 0.13 + j], [low, high], color=colors[dataset],
                         linewidth=0.55, alpha=0.35)
        axes[2].plot([xpos - 0.13, xpos + 0.13],
                     [sub.lower_coverage_mean_effect_pp.mean(), sub.upper_coverage_mean_effect_pp.mean()],
                     color=colors[dataset], linewidth=2.0, solid_capstyle="round")
    axes[2].axhline(0, color=OKABE["black"], linewidth=0.7)
    axes[2].set_xticks([0, 1], ["GSE294385\n(n=5)", "Valdeolivas\n(n=6)"])
    axes[2].set_ylabel("Effect (percentage points)")
    axes[2].set_title("Coverage halves", loc="left", fontweight="bold")
    axes[2].legend(frameon=False, fontsize=6.5, loc="upper center", ncol=2)
    axes[2].text(0.02, -0.28, "Dependent accepted shifts; no inferential P value",
                 transform=axes[2].transAxes, fontsize=6.3, color=OKABE["gray"])

    # Each pair appears twice (IS-low/high) in S14b: count an accepted offset once.
    keys = ["dataset", "sample", "patient_id", "shift_array_row", "shift_array_col"]
    passing_sections = sections.loc[sections.section_status.eq("PASS"), ["dataset", "sample"]]
    shift_rows = shift_rows.merge(passing_sections, on=["dataset", "sample"], how="inner", validate="many_to_one")
    pair_rows = pair_rows.merge(passing_sections, on=["dataset", "sample"], how="inner", validate="many_to_one")
    passed = pair_rows.loc[pair_rows.status.eq("PASS") & pair_rows.pair_passed_minimum].copy()
    passed = passed.drop_duplicates(keys + ["pathology_pair"])
    accepted = shift_rows.loc[shift_rows.status.eq("PASS")].groupby("dataset").size()
    assert int(accepted["GSE294385"]) == 672 and int(accepted["Valdeolivas_2024"]) == 723
    pair_order = ["tumor->tumor", "tumor->stroma", "stroma->tumor", "stroma->stroma"]
    patient_support = []
    accepted_patient = shift_rows.loc[shift_rows.status.eq("PASS")].groupby(["dataset", "patient_id"])
    for (dataset, patient_id), accepted_rows in accepted_patient:
        for pair in pair_order:
            subset = passed.loc[passed.dataset.eq(dataset) & passed.patient_id.eq(patient_id) &
                                passed.pathology_pair.eq(pair)]
            patient_support.append({"dataset": dataset, "patient_id": patient_id,
                                    "pathology_pair": pair,
                                    "accepted_shift_section_n": len(accepted_rows),
                                    "qualified_shift_section_n": len(subset),
                                    "qualified_section_n": int(subset["sample"].nunique()),
                                    "qualified_fraction": len(subset) / len(accepted_rows)})
    patient_support_df = pd.DataFrame(patient_support)
    patient_support_df.to_csv(audit_dir / "Table_S14g_patient_pair_support.csv", index=False)
    support = []
    for dataset in ["GSE294385", "Valdeolivas_2024"]:
        for pair in pair_order:
            subset = passed.loc[passed.dataset.eq(dataset) & passed.pathology_pair.eq(pair)]
            support.append({"dataset": dataset, "pathology_pair": pair,
                            "accepted_shift_section_n": int(accepted[dataset]),
                            "qualified_shift_section_n": len(subset),
                            "qualified_patient_n": int(subset.patient_id.nunique()),
                            "qualified_fraction": len(subset) / int(accepted[dataset])})
    support_df = pd.DataFrame(support)
    support_df.to_csv(audit_dir / "Table_S14f_pathology_pair_support.csv", index=False)
    assert (support_df.loc[(support_df.dataset == "Valdeolivas_2024") &
                           support_df.pathology_pair.isin(["tumor->stroma", "stroma->tumor"]),
                           "qualified_shift_section_n"] == 0).all()
    x_pair = np.arange(4)
    for offset, dataset in zip((-width / 2, width / 2), ["GSE294385", "Valdeolivas_2024"]):
        sub = support_df.loc[support_df.dataset.eq(dataset)].set_index("pathology_pair").loc[pair_order]
        axes[3].bar(x_pair + offset, 100 * sub.qualified_fraction.to_numpy(), width,
                    color=colors[dataset], edgecolor="white", linewidth=0.35)
    axes[3].set_xticks(x_pair, ["T→T", "T→S", "S→T", "S→S"])
    axes[3].tick_params(axis="x", labelsize=8)
    axes[3].set_ylabel("Accepted offsets with pair (%)")
    axes[3].set_ylim(0, 100)
    axes[3].set_title("Pathology-pair support", loc="left", fontweight="bold")
    # Zero qualifying offsets are technical support failures, not zero-valued
    # cross-compartment effects. Label the two external empty positions in-plot.
    for pair_index in (1, 2):
        axes[3].text(x_pair[pair_index] + width / 2, 4.0, "0/723\nNE",
                     ha="center", va="bottom", fontsize=8,
                     fontweight="bold", color=colors["Valdeolivas_2024"])
    axes[3].text(0.02, -0.28, "Technical offsets only; shifts are not patients",
                 transform=axes[3].transAxes, fontsize=7, color=OKABE["gray"])

    for axis, label in zip(axes, "ABCD"):
        axis.text(-0.10, 1.08, label, transform=axis.transAxes, fontsize=10, fontweight="bold", va="top")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(False)
    fig.tight_layout(pad=0.7)
    for ext in ("pdf", "png"):
        fig.savefig(figure_dir / f"Fig3_v11_10_edge_selection_audit.{ext}", dpi=610 if ext == "png" else None,
                    facecolor="white")
    manifest = {"figure": "Figure 3", "inputs": {
        "cohort": sha256(audit_dir / "Table_patient_equal_cohort_selection.csv"),
        "patients": sha256(audit_dir / "Table_patient_selection_summary.csv"),
        "halves": sha256(audit_dir / "Table_patient_coverage_halves.csv"),
        "shift_rows": sha256(audit_dir / "Table_shift_edge_selection.csv"),
        "pair_rows": sha256(audit_dir / "Table_shift_pair_and_edge_hops.csv"),
        "passing_sections": sha256(audit_dir / "Table_section_selection_summary.csv")},
        "script": sha256(Path(__file__)), "n_patients": len(patients),
        "note": "Descriptive post hoc diagnostic; no shift is an independent unit and no P value is shown."}
    (figure_dir / "Fig3_v11_10_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(figure_dir / "Fig3_v11_10_edge_selection_audit.pdf")


if __name__ == "__main__":
    main()
