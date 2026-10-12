#!/usr/bin/env python
"""Analysis: patient-aware spatial-shift and same-spot-confounding audit.

Date: 2026-10-05
Random seed: 42
Python and key package versions are recorded in the output manifest.

The shift analysis translates the *intact* SS label field on the Visium array
lattice. It is a spatially structured negative control, not an exact permutation
test: irregular tissue boundaries and Layer3 matching change the eligible mask.
"""

from __future__ import annotations

import argparse
import io
import json
import tarfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree


np.random.seed(42)


def lattice_positions(archive: Path) -> pd.DataFrame:
    with tarfile.open(archive, "r:gz") as handle:
        members = [member for member in handle if member.name.endswith("/spatial/tissue_positions.csv")]
        if len(members) != 1:
            raise ValueError(f"Expected one tissue_positions.csv in {archive}")
        stream = handle.extractfile(members[0])
        if stream is None:
            raise ValueError(f"Unreadable positions file in {archive}")
        positions = pd.read_csv(io.BytesIO(stream.read()))
    needed = {"barcode", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"}
    if not needed.issubset(positions):
        raise ValueError(f"Position columns missing in {archive}")
    if positions.barcode.duplicated().any():
        raise ValueError(f"Duplicate barcodes in {archive}")
    return positions[list(needed)]


def select_offsets() -> list[tuple[int, int, int]]:
    """Sixteen fixed-seed translations at each hex-lattice radius 5 through 10."""
    rng = np.random.default_rng(42)
    selected = []
    for radius in range(5, 11):
        choices = []
        for dr in range(-radius, radius + 1):
            for dc in range(-2 * radius, 2 * radius + 1):
                if (dr + dc) % 2:
                    continue
                hex_distance = max(abs(dr), abs((dc + dr) // 2), abs((dc - dr) // 2))
                if hex_distance == radius:
                    choices.append((dr, dc, radius))
        indices = rng.choice(len(choices), size=16, replace=False)
        selected.extend(choices[index] for index in indices)
    if len(set(selected)) != 96:
        raise AssertionError("Spatial translations not unique")
    return selected


def matched_layer_contrast(
    source: np.ndarray, target: np.ndarray, immune: np.ndarray,
    stroma: np.ndarray, layer: np.ndarray,
) -> tuple[float, float, int, int]:
    """IS-high minus IS-low near-SS prevalence among SS-low source spots.

    Compare within the same (source Layer3, target Layer3) pair, weighting each
    pair by the smaller of its two group edge counts. No edge is an inference unit.
    """
    eligible = ~stroma[source]
    table = pd.DataFrame({
        "source_layer": layer[source[eligible]],
        "target_layer": layer[target[eligible]],
        "IS_high": immune[source[eligible]],
        "neighbor_SS_high": stroma[target[eligible]].astype(np.int8),
    })
    grouped = table.groupby(["source_layer", "target_layer", "IS_high"], observed=True)[
        "neighbor_SS_high"
    ].agg(["sum", "count"])
    cells = grouped.unstack("IS_high", fill_value=0)
    if not {("count", False), ("count", True)}.issubset(cells.columns):
        return np.nan, 0.0, 0, 0
    sufficient = (cells[("count", False)] >= 30) & (cells[("count", True)] >= 30)
    cells = cells.loc[sufficient]
    if cells.empty:
        return np.nan, 0.0, 0, 0
    high_rate = cells[("sum", True)] / cells[("count", True)]
    low_rate = cells[("sum", False)] / cells[("count", False)]
    weights = np.minimum(cells[("count", True)], cells[("count", False)])
    contrast = float(np.average(high_rate - low_rate, weights=weights))
    total_high_edges = int((eligible & immune[source]).sum())
    matched_high_edges = int(cells[("count", True)].sum())
    return contrast, matched_high_edges / total_high_edges, len(cells), matched_high_edges


def shifted_matched_layer_contrast(
    source: np.ndarray, target: np.ndarray, mapped: np.ndarray,
    compatible: np.ndarray, immune: np.ndarray, stroma: np.ndarray,
    layer: np.ndarray,
) -> tuple[float, float, int]:
    """Observed-minus-shifted IS contrast, restricted to SS-low source spots.

    Both contrasts use the identical compatible edges and source/target Layer3
    pairs, so the result addresses same-spot SS status and broad SS field
    structure simultaneously. It remains a sensitivity contrast, not a causal
    or exact-randomization test.
    """
    eligible = compatible & ~stroma[source]
    if int(eligible.sum()) < 100:
        return np.nan, 0.0, 0
    table = pd.DataFrame({
        "source_layer": layer[source[eligible]],
        "target_layer": layer[target[eligible]],
        "IS_high": immune[source[eligible]],
        "observed_neighbor": stroma[target[eligible]].astype(np.int8),
        "shifted_neighbor": stroma[mapped[eligible]].astype(np.int8),
    })
    grouped = table.groupby(["source_layer", "target_layer", "IS_high"], observed=True).agg(
        n=("observed_neighbor", "size"),
        observed_sum=("observed_neighbor", "sum"),
        shifted_sum=("shifted_neighbor", "sum"),
    )
    cells = grouped.unstack("IS_high", fill_value=0)
    if not {("n", False), ("n", True)}.issubset(cells.columns):
        return np.nan, 0.0, 0
    cells = cells.loc[(cells[("n", False)] >= 30) & (cells[("n", True)] >= 30)]
    if cells.empty:
        return np.nan, 0.0, 0
    observed_contrast = cells[("observed_sum", True)] / cells[("n", True)] - (
        cells[("observed_sum", False)] / cells[("n", False)]
    )
    shifted_contrast = cells[("shifted_sum", True)] / cells[("n", True)] - (
        cells[("shifted_sum", False)] / cells[("n", False)]
    )
    weights = np.minimum(cells[("n", True)], cells[("n", False)])
    all_high_SS_low_edges = int((~stroma[source] & immune[source]).sum())
    matched_high = int(cells[("n", True)].sum())
    return float(np.average(observed_contrast - shifted_contrast, weights=weights)), matched_high / all_high_SS_low_edges, len(cells)


def one_section(record: pd.Series, root: Path, offsets: list[tuple[int, int, int]]) -> tuple[dict, list[dict]]:
    gsm = record["sample"]
    scores = pd.read_csv(root / "sections" / f"Table_{gsm}_spot_scores_and_boundaries.csv")
    pathology = pd.read_csv(
        root / "annotations" / f"Table_{gsm}_author_pathology_spot_compartments.csv",
        usecols=["barcode", "Layer3"],
    )
    positions = lattice_positions(Path(record["archive"]))
    data = scores.merge(pathology, on="barcode", how="left", validate="one_to_one")
    data = data.merge(positions, on="barcode", how="left", validate="one_to_one")
    if len(data) != int(record["n_spots"]) or data[
        ["Layer3", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    ].isna().any().any():
        raise ValueError(f"Spot, pathology, or lattice mapping failure for {gsm}")
    if not np.array_equal(data.row.to_numpy(int), data.pxl_row_in_fullres.to_numpy(int)) or not np.array_equal(
        data.col.to_numpy(int), data.pxl_col_in_fullres.to_numpy(int)
    ):
        raise ValueError(f"Pixel coordinate disagreement for {gsm}")
    array_row = data.array_row.to_numpy(int)
    array_col = data.array_col.to_numpy(int)
    if np.any((array_row + array_col) % 2):
        raise ValueError(f"Non-Visium lattice parity for {gsm}")
    layer = data.Layer3.to_numpy(str)
    immune = data.IS_high_q75.to_numpy(bool)
    stroma = data.SS_high_q75.to_numpy(bool)
    xy = data[["row", "col"]].to_numpy(float)
    distance, neighbor = cKDTree(xy).query(xy, k=7)
    pitch = float(np.median(distance[:, 1]))
    valid = (distance[:, 1:] <= 1.6 * pitch).reshape(-1)
    source = np.repeat(np.arange(len(data)), 6)[valid]
    target = neighbor[:, 1:].reshape(-1)[valid]
    if len(source) < 100:
        raise ValueError(f"Too few local graph edges for {gsm}")
    conditional, conditional_coverage, n_strata, matched_edges = matched_layer_contrast(
        source, target, immune, stroma, layer
    )
    selected = immune[source]
    selected_target = target[selected]
    if len(selected_target) < 100:
        raise ValueError(f"Too few local IS-high edges for {gsm}")

    lookup = np.full((int(array_row.max()) + 1, int(array_col.max()) + 1), -1, dtype=np.int32)
    lookup[array_row, array_col] = np.arange(len(data), dtype=np.int32)
    if int((lookup >= 0).sum()) != len(data):
        raise ValueError(f"Duplicate array coordinates for {gsm}")
    shift_rows = []
    for dr, dc, radius in offsets:
        rr = array_row[target] + dr
        cc = array_col[target] + dc
        in_grid = (rr >= 0) & (rr < lookup.shape[0]) & (cc >= 0) & (cc < lookup.shape[1])
        mapped = np.full(len(target), -1, dtype=np.int32)
        mapped[in_grid] = lookup[rr[in_grid], cc[in_grid]]
        compatible = mapped >= 0
        compatible[compatible] &= layer[mapped[compatible]] == layer[target[compatible]]
        high_compatible = compatible & selected
        n_matched = int(high_compatible.sum())
        coverage = n_matched / len(selected_target)
        if n_matched < 100 or coverage < 0.10:
            continue
        observed = float(stroma[target[high_compatible]].mean())
        shifted = float(stroma[mapped[high_compatible]].mean())
        conditional_did, did_coverage, did_pairs = shifted_matched_layer_contrast(
            source, target, mapped, compatible, immune, stroma, layer
        )
        shift_rows.append({
            "sample": gsm, "patient_id": record["patient_id"], "tissue_context": record["tissue_context"],
            "shift_array_row": dr, "shift_array_col": dc, "hex_radius": radius,
            "n_matched_directed_edges": n_matched, "matched_edge_fraction": coverage,
            "observed_SS_high_fraction_on_matched_edges": observed,
            "shifted_SS_high_fraction_on_matched_edges": shifted,
            "observed_minus_shifted_fraction": observed - shifted,
            "SS_low_source_shifted_matched_difference_in_differences": conditional_did,
            "SS_low_source_shifted_matched_IS_high_edge_fraction_of_original": did_coverage,
            "SS_low_source_shifted_matched_Layer3_pairs": did_pairs,
        })
    if len(shift_rows) < 24:
        raise ValueError(f"Only {len(shift_rows)} valid spatial shifts for {gsm}")
    shifts = pd.DataFrame(shift_rows)
    if shifts.SS_low_source_shifted_matched_difference_in_differences.notna().sum() < 24:
        raise ValueError(f"Too few valid conditional shifts for {gsm}")
    summary = {
        "sample": gsm, "patient_id": record["patient_id"], "tissue_context": record["tissue_context"],
        "n_spots": len(data), "n_short_directed_edges": len(source),
        "n_IS_high_short_directed_edges": len(selected_target),
        "same_spot_SS_high_fraction_among_IS_high": float(stroma[immune].mean()),
        "n_fixed_candidate_shifts": len(offsets), "n_accepted_shifts": len(shifts),
        "shift_matched_edge_fraction_median": float(shifts.matched_edge_fraction.median()),
        "shift_control_delta_mean": float(shifts.observed_minus_shifted_fraction.mean()),
        "shift_control_delta_median": float(shifts.observed_minus_shifted_fraction.median()),
        "fraction_shifts_observed_greater_than_shifted": float((shifts.observed_minus_shifted_fraction > 0).mean()),
        "SS_low_source_shifted_matched_difference_in_differences_mean": float(
            shifts.SS_low_source_shifted_matched_difference_in_differences.mean()
        ),
        "SS_low_source_shifted_matched_IS_high_edge_fraction_of_original_median": float(
            shifts.SS_low_source_shifted_matched_IS_high_edge_fraction_of_original.median()
        ),
        "SS_low_source_IS_high_minus_IS_low_neighbor_SS_fraction": conditional,
        "SS_low_source_matched_IS_high_edge_fraction": conditional_coverage,
        "SS_low_source_n_matched_Layer3_pairs": n_strata,
        "SS_low_source_n_matched_IS_high_edges": matched_edges,
    }
    return summary, shift_rows


def patient_summary(frame: pd.DataFrame, metric: str, scenario: str, rng: np.random.Generator) -> tuple[dict, pd.DataFrame]:
    patients = frame.groupby("patient_id", as_index=False)[metric].mean().dropna()
    values = patients[metric].to_numpy(float)
    if len(values) < 2:
        raise ValueError(f"Too few patients for {scenario}/{metric}")
    replicates = rng.choice(values, size=(10000, len(values)), replace=True).mean(axis=1)
    lower, upper = np.quantile(replicates, [0.025, 0.975])
    row = {
        "scenario": scenario, "metric": metric, "n_sections": len(frame), "n_patients": len(values),
        "n_positive_sections": int((frame[metric] > 0).sum()),
        "n_positive_patients": int((values > 0).sum()),
        "estimate": float(values.mean()), "ci_lower": float(lower), "ci_upper": float(upper),
        "ci_method": "patient bootstrap percentile, 10000 draws, seed 42",
        "inference_status": "exploratory; spatial shifts are structured negative controls, not exact exchangeable permutations",
    }
    patients["scenario"] = scenario
    return row, patients


def draw_patient_figure(patients: pd.DataFrame, summary: pd.DataFrame, out: Path, style: Path | None) -> None:
    if style and style.is_file():
        plt.style.use(style)
    all_patients = patients.loc[patients.scenario.eq("all")]
    metrics = [
        ("shift_control_delta_mean", "Observed minus shifted SS field"),
        ("SS_low_source_IS_high_minus_IS_low_neighbor_SS_fraction", "IS-high minus IS-low, SS-low sources"),
        ("SS_low_source_shifted_matched_difference_in_differences_mean", "Both controls: observed minus shifted"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.5), sharey=True)
    order = sorted(all_patients.patient_id.unique())
    for axis, (metric, title) in zip(axes, metrics):
        data = all_patients.loc[all_patients.metric.eq(metric)].set_index("patient_id").reindex(order)
        axis.scatter(data[metric].to_numpy(float) * 100, np.arange(len(order)), color="#0072B2", zorder=3)
        pooled = summary.loc[summary.scenario.eq("all") & summary.metric.eq(metric)].iloc[0]
        axis.axvline(0, color="grey", linestyle="--", linewidth=1)
        axis.axvspan(pooled.ci_lower * 100, pooled.ci_upper * 100, color="#D55E00", alpha=0.15)
        axis.axvline(pooled.estimate * 100, color="#D55E00", linewidth=1.5)
        axis.set_title(title)
        axis.set_xlabel("Difference in SS-high neighbours (percentage points)")
        axis.set_yticks(np.arange(len(order)), order)
        axis.invert_yaxis()
    fig.suptitle("GSE294385: patient-level spatial sensitivity")
    fig.savefig(out / "Fig_GSE294385_spatial_shift_and_SS_low_source_audit.pdf")
    fig.savefig(out / "Fig_GSE294385_spatial_shift_and_SS_low_source_audit.png", dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--style", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    sections = pd.read_csv(args.root / "Table_GSE294385_37_section_effects.csv")
    if len(sections) != 37 or sections.patient_id.nunique() != 11:
        raise ValueError("Expected 37 sections from 11 patients")
    offsets = select_offsets()
    rows = []
    all_shifts = []
    for _, record in sections.iterrows():
        one, shifts = one_section(record, args.root, offsets)
        rows.append(one)
        all_shifts.extend(shifts)
        print(f"{one['sample']}: {one['n_accepted_shifts']} accepted shifts; "
              f"shift delta {one['shift_control_delta_mean']:.4f}; "
              f"SS-low contrast {one['SS_low_source_IS_high_minus_IS_low_neighbor_SS_fraction']:.4f}", flush=True)
    section_table = pd.DataFrame(rows)
    section_table.to_csv(args.out / "Table_GSE294385_spatial_autocorr_audit_by_section.csv", index=False)
    pd.DataFrame(all_shifts).to_csv(args.out / "Table_GSE294385_spatial_shift_diagnostics.csv", index=False)
    rng = np.random.default_rng(42)
    summaries = []
    patient_frames = []
    groups = [("all", section_table)] + [(key, group) for key, group in section_table.groupby("tissue_context")]
    for scenario, frame in groups:
        for metric in (
            "shift_control_delta_mean",
            "SS_low_source_IS_high_minus_IS_low_neighbor_SS_fraction",
            "SS_low_source_shifted_matched_difference_in_differences_mean",
        ):
            row, patients = patient_summary(frame, metric, scenario, rng)
            summaries.append(row)
            patients["metric"] = metric
            patient_frames.append(patients)
    summary = pd.DataFrame(summaries)
    patient_values = pd.concat(patient_frames, ignore_index=True)
    summary.to_csv(args.out / "Table_GSE294385_spatial_autocorr_audit_patient_summary.csv", index=False)
    patient_values.to_csv(args.out / "Table_GSE294385_spatial_autocorr_audit_patient_values.csv", index=False)
    draw_patient_figure(patient_values, summary, args.out, args.style)
    manifest = {
        "date": "2026-10-05", "analysis": "spatial translation and SS-low source matching",
        "seed": 42, "python": __import__("sys").version, "numpy": np.__version__,
        "pandas": pd.__version__, "scipy": scipy.__version__,
        "n_sections": len(sections), "n_patients": int(sections.patient_id.nunique()),
        "translation": "96 fixed Visium lattice offsets, 16 each at hex radius 5-10; intact SS field translated; same target Layer3 and >=10% edge overlap",
        "conditional_contrast": "IS-high vs IS-low among SS-low source spots, same source and target Layer3 pair, >=30 edges/group/pair, min-count weighting",
        "joint_control": "On identically eligible edges and Layer3 pairs, (observed IS-high minus IS-low) minus (shifted-field IS-high minus IS-low), SS-low source spots only",
        "inference_unit": "patient, after section means",
        "limitation": "Not exact randomization: offsets correlate, tissue/pathology masks vary; no causal or direct cell-contact inference",
    }
    (args.out / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (args.out / "_analysis_outputs.md").write_text(
        "# Analysis outputs\nGenerated: 2026-10-05\nStudy type: exploratory patient-aware Visium spatial sensitivity\n\n"
        "## Tables\n- Table_GSE294385_spatial_autocorr_audit_by_section.csv\n"
        "- Table_GSE294385_spatial_shift_diagnostics.csv\n"
        "- Table_GSE294385_spatial_autocorr_audit_patient_summary.csv\n"
        "- Table_GSE294385_spatial_autocorr_audit_patient_values.csv\n\n"
        "## Figures\n- Fig_GSE294385_spatial_shift_and_SS_low_source_audit.pdf / .png\n",
        encoding="utf-8",
    )
    print("\nPatient-level summary (fractions; multiply by 100 for percentage points):")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
