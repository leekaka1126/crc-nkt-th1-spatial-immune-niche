#!/usr/bin/env python
"""Patient-level effect plot for the locked CRC spatial analysis.

Key message: primary-site patient effects vary and both confidence
intervals include zero; tissue spots and replicate sections are not plotted as
independent observations. Intended for 30–120-second journal reading.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {"Valdeolivas_2024": "#0072B2", "GSE294385": "#D55E00"}


def draw_panel(ax, patients: pd.DataFrame, summary: pd.Series, label: str, color: str) -> None:
    patients = patients.sort_values("patient_id").reset_index(drop=True)
    effects = patients.pathology_matched_shift_effect_fraction.to_numpy(float) * 100
    y = np.arange(len(patients), 0, -1, dtype=float)
    ax.scatter(effects, y, s=40, marker="o", facecolor=color, edgecolor="black", linewidth=0.45, zorder=3)
    estimate = float(summary.estimate_fraction) * 100
    lower = float(summary.ci_lower_fraction) * 100
    upper = float(summary.ci_upper_fraction) * 100
    ax.hlines(0, lower, upper, color="black", linewidth=1.5, zorder=2)
    ax.plot([lower, lower], [-0.10, 0.10], color="black", linewidth=1)
    ax.plot([upper, upper], [-0.10, 0.10], color="black", linewidth=1)
    ax.scatter([estimate], [0], s=62, marker="D", facecolor="black", edgecolor="black", zorder=4)
    ax.axvline(0, color="#666666", linewidth=0.9, linestyle=(0, (3, 3)), zorder=1)
    ax.set_yticks([*y, 0], [*patients.patient_id.tolist(), "Mean"])
    ax.set_ylim(-0.65, len(patients) + 0.8)
    ax.set_xlim(-5.5, 5.5)
    ax.set_xticks([-4, -2, 0, 2, 4])
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", alpha=0.25)
    ax.grid(axis="y", visible=False)
    ax.set_title(label, loc="left", fontsize=10, pad=6)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    patients = pd.read_csv(args.run_dir / "Table_patient_effects.csv", dtype={"patient_id": str})
    summary = pd.read_csv(args.run_dir / "Table_patient_equal_summary.csv")
    scenarios = [
        ("Valdeolivas_2024", "all", "A  External primary CRC (6 patients)"),
        ("GSE294385", "Primary_CRC", "B  GSE294385 primary CRC (5 patients)"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.8), sharex=True)
    for ax, (dataset, context, label) in zip(axes, scenarios):
        subset = patients.loc[patients.dataset.eq(dataset) & patients.context.eq(context)]
        result = summary.loc[summary.dataset.eq(dataset) & summary.context.eq(context)]
        if len(result) != 1 or len(subset) != int(result.iloc[0].n_eligible_patients):
            raise ValueError(f"Patient denominator mismatch for {dataset}/{context}")
        draw_panel(ax, subset, result.iloc[0], label, COLORS[dataset])
    fig.supxlabel("Pathology-matched contrast (percentage points)", y=0.035, fontsize=9)
    fig.subplots_adjust(left=0.13, right=0.99, bottom=0.20, top=0.89, wspace=0.43)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=610, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
