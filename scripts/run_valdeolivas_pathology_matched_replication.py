#!/usr/bin/env python
"""Frozen-panel, pathology-matched CRC spatial replication with patient-level summaries.

Run only after the effect-blind category/QC lock. All attempted sections and all
96 fixed translations for each section reaching that stage are retained. A translation is a structured
negative control, not an exchangeable permutation or an independent sample.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import zipfile

import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
PANEL_SCRIPT = ROOT / "manuscript/CRC_NKT_v11_1/code_release/scripts/05_spatial_neighborhood_analysis.py"
SHIFT_SCRIPT = ROOT / "scripts/gse294385_spatial_autocorr_audit.py"
JOIN_SCRIPT = ROOT / "scripts/audit_valdeolivas_spot_joins.py"
STRICT = {"tumor", "stroma"}


def import_file(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


score_mod = import_file("valdeolivas_fixed_score", PANEL_SCRIPT)
shift_mod = import_file("valdeolivas_fixed_shift", SHIFT_SCRIPT)
join_mod = import_file("valdeolivas_input_join", JOIN_SCRIPT)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def true_flags(series: pd.Series, label: str) -> np.ndarray:
    cleaned = series.astype(str).str.lower()
    if not cleaned.isin(["true", "false"]).all():
        raise ValueError(f"Missing or nonboolean {label}")
    return cleaned.eq("true").to_numpy(bool)


def matched_effect(data: pd.DataFrame, offsets: list[tuple[int, int, int]]) -> tuple[dict, list[dict]]:
    """Apply original same-pair and intact-field shift logic to strict T/S edges."""
    needed = ["array_row", "array_col", "row", "col", "compartment", "IS_high_q75", "SS_high_q75"]
    if data[needed].isna().any().any():
        raise ValueError("Missing coordinate, pathology, or fixed-score field")
    ar, ac = data.array_row.to_numpy(int), data.array_col.to_numpy(int)
    if ((ar + ac) % 2).any() or data[["array_row", "array_col"]].duplicated().any():
        raise ValueError("Invalid or repeated Visium lattice address")
    layer = data.compartment.astype(str).to_numpy()
    if set(layer) - {"tumor", "stroma", "other"}:
        raise ValueError("Unmapped pathology category")
    immune = true_flags(data.IS_high_q75, "IS_high_q75")
    stroma = true_flags(data.SS_high_q75, "SS_high_q75")
    xy = data[["row", "col"]].to_numpy(float)
    if not np.isfinite(xy).all() or len(data) < 7:
        raise ValueError("Insufficient finite spot coordinates")
    distance, neighbor = cKDTree(xy).query(xy, k=7)
    pitch = float(np.median(distance[:, 1]))
    valid = (distance[:, 1:] <= 1.6 * pitch).reshape(-1)
    all_source = np.repeat(np.arange(len(data)), 6)[valid]
    all_target = neighbor[:, 1:].reshape(-1)[valid]
    strict_mask = np.isin(layer[all_source], list(STRICT)) & np.isin(layer[all_target], list(STRICT))
    source, target = all_source[strict_mask], all_target[strict_mask]
    pair_counts = Counter(f"{layer[i]}->{layer[j]}" for i, j in zip(source, target))
    selected = immune[source]
    n_selected = int(selected.sum())
    base = {
        "n_spots": len(data),
        "n_all_short_directed_edges": len(all_source),
        "n_strict_short_directed_edges": len(source),
        "n_tumor_to_stroma_directed_edges": pair_counts.get("tumor->stroma", 0),
        "n_IS_high_strict_edges": n_selected,
        "median_nearest_neighbor_pitch_fullres_px": pitch,
        "same_spot_SS_high_fraction_among_IS_high": float(stroma[immune].mean()) if immune.any() else np.nan,
        "n_fixed_candidate_shifts": len(offsets),
    }
    if len(source) < 100 or n_selected < 100:
        return {**base, "status": "FAIL", "reason": "fewer than 100 strict local or IS-high strict edges"}, []
    n_ss_low_high = int((~stroma[source] & immune[source]).sum())
    n_ss_low_low = int((~stroma[source] & ~immune[source]).sum())
    base["n_SS_low_IS_high_strict_edges"] = n_ss_low_high
    base["n_SS_low_IS_low_strict_edges"] = n_ss_low_low
    if min(n_ss_low_high, n_ss_low_low) < 30:
        return {**base, "status": "FAIL", "reason": "fewer than 30 SS-low strict edges in one IS group"}, []
    conditional, initial_coverage, initial_pairs, initial_high = shift_mod.matched_layer_contrast(
        source, target, immune, stroma, layer
    )
    base.update({
        "initial_matched_IS_contrast": conditional,
        "initial_matched_IS_high_edge_coverage": initial_coverage,
        "initial_matched_pathology_pairs": initial_pairs,
        "initial_matched_IS_high_edges": initial_high,
    })
    lookup = np.full((int(ar.max()) + 1, int(ac.max()) + 1), -1, dtype=np.int32)
    lookup[ar, ac] = np.arange(len(data), dtype=np.int32)
    if int((lookup >= 0).sum()) != len(data):
        raise ValueError("Lattice lookup loses spots")
    shift_rows = []
    for dr, dc, radius in offsets:
        rr, cc = ar[target] + dr, ac[target] + dc
        in_grid = (rr >= 0) & (rr < lookup.shape[0]) & (cc >= 0) & (cc < lookup.shape[1])
        mapped = np.full(len(target), -1, dtype=np.int32)
        mapped[in_grid] = lookup[rr[in_grid], cc[in_grid]]
        compatible = mapped >= 0
        compatible[compatible] &= layer[mapped[compatible]] == layer[target[compatible]]
        n_matched = int((compatible & selected).sum())
        coverage = n_matched / n_selected
        record = {
            "shift_array_row": dr, "shift_array_col": dc, "hex_radius": radius,
            "n_matched_IS_high_edges_before_pair_filter": n_matched,
            "matched_IS_high_edge_fraction_before_pair_filter": coverage,
            "n_matched_pathology_pairs": 0,
            "matched_IS_high_edge_fraction_after_pair_filter": 0.0,
            "observed_minus_shifted_IS_contrast": np.nan,
        }
        if n_matched < 100 or coverage < 0.10:
            record["status"] = "FAIL_COVERAGE"
        else:
            effect, pair_coverage, n_pairs = shift_mod.shifted_matched_layer_contrast(
                source, target, mapped, compatible, immune, stroma, layer
            )
            record.update({
                "n_matched_pathology_pairs": n_pairs,
                "matched_IS_high_edge_fraction_after_pair_filter": pair_coverage,
                "observed_minus_shifted_IS_contrast": effect,
            })
            record["status"] = "PASS" if np.isfinite(effect) and n_pairs >= 1 else "FAIL_PAIR_MINIMUM"
        shift_rows.append(record)
    valid_shifts = [record for record in shift_rows if record["status"] == "PASS"]
    base.update({
        "n_accepted_shifts": len(valid_shifts),
        "n_coverage_fail_shifts": sum(record["status"] == "FAIL_COVERAGE" for record in shift_rows),
        "n_pair_minimum_fail_shifts": sum(record["status"] == "FAIL_PAIR_MINIMUM" for record in shift_rows),
        "matched_IS_high_edge_fraction_after_pair_filter_median": float(np.median([
            record["matched_IS_high_edge_fraction_after_pair_filter"] for record in valid_shifts
        ])) if valid_shifts else np.nan,
    })
    if len(valid_shifts) < 24:
        return {**base, "status": "FAIL", "reason": f"only {len(valid_shifts)} of 96 translations passed fixed thresholds"}, shift_rows
    return {
        **base,
        "status": "PASS",
        "reason": "",
        "pathology_matched_shift_effect_fraction": float(np.mean([
            record["observed_minus_shifted_IS_contrast"] for record in valid_shifts
        ])),
    }, shift_rows


def external_data(section_path: Path, annotation_zip: zipfile.ZipFile, mapping: dict[str, str], scratch: Path) -> tuple[pd.DataFrame, dict]:
    sample = join_mod.SAMPLE_RE.search(section_path.name)
    if not sample:
        raise ValueError(f"Unparseable sample file: {section_path.name}")
    patient_id, rep = sample.group(1), int(sample.group(2))
    _, labels = join_mod.annotation_for_sample(annotation_zip, patient_id, rep)
    with zipfile.ZipFile(section_path) as archive:
        h5_member = join_mod.exactly_one(archive, "/filtered_feature_bc_matrix.h5")
        with tempfile.TemporaryDirectory(prefix="visium_", dir=scratch) as temp_name:
            h5_path = Path(temp_name) / "filtered_feature_bc_matrix.h5"
            with archive.open(h5_member) as source, h5_path.open("wb") as target:
                shutil.copyfileobj(source, target, length=1024 * 1024)
            expression = score_mod.read_10x_panel(h5_path, score_mod.IMMUNE_GENES + score_mod.MATRIX_GENES)
        _, positions = join_mod.position_metadata(archive)
    tissue_positions = {barcode: value for barcode, value in positions.items() if value[0] == 1}
    if set(expression.index) != set(labels) or set(expression.index) != set(tissue_positions):
        raise ValueError("Expression, annotation and in-tissue coordinate barcodes do not match")
    data = pd.DataFrame(index=expression.index)
    data.index.name = "barcode"
    data["array_row"] = [tissue_positions[barcode][1] for barcode in data.index]
    data["array_col"] = [tissue_positions[barcode][2] for barcode in data.index]
    data["row"] = [tissue_positions[barcode][3] for barcode in data.index]
    data["col"] = [tissue_positions[barcode][4] for barcode in data.index]
    data["compartment"] = [mapping[labels[barcode]] for barcode in data.index]
    data["IS_within_section"] = score_mod.zmean(expression, score_mod.IMMUNE_GENES)
    data["SS_within_section"] = score_mod.zmean(expression, score_mod.MATRIX_GENES)
    data["IS_high_q75"] = data.IS_within_section >= data.IS_within_section.quantile(0.75)
    data["SS_high_q75"] = data.SS_within_section >= data.SS_within_section.quantile(0.75)
    variable = {
        "n_variable_IS_genes": int((expression[score_mod.IMMUNE_GENES].std(axis=0, ddof=1) > 0).sum()),
        "n_variable_SS_genes": int((expression[score_mod.MATRIX_GENES].std(axis=0, ddof=1) > 0).sum()),
        "IS_high_spot_fraction": float(data.IS_high_q75.mean()),
        "SS_high_spot_fraction": float(data.SS_high_q75.mean()),
        "IS_score_sd": float(data.IS_within_section.std(ddof=1)),
        "SS_score_sd": float(data.SS_within_section.std(ddof=1)),
    }
    return data, variable


def source_data(record: pd.Series, source_spots: Path, source_annotations: Path) -> pd.DataFrame:
    gsm = str(record["sample"])
    scores = pd.read_csv(source_spots / f"Table_{gsm}_spot_scores_and_boundaries.csv")
    labels = pd.read_csv(source_annotations / f"Table_{gsm}_author_pathology_spot_compartments.csv",
                         usecols=["barcode", "compartment"])
    positions = shift_mod.lattice_positions(Path(record["archive"]))
    if scores.barcode.duplicated().any() or labels.barcode.duplicated().any() or positions.barcode.duplicated().any():
        raise ValueError("Duplicate source barcode")
    source_label = scores[["barcode", "compartment"]].rename(columns={"compartment": "score_compartment"})
    source_join = source_label.merge(labels, on="barcode", how="outer", validate="one_to_one", indicator=True)
    if not source_join._merge.eq("both").all() or not source_join.score_compartment.eq(source_join.compartment).all():
        raise ValueError("Precomputed score and author pathology compartments disagree")
    scores = scores.drop(columns=["compartment"])
    data = scores.merge(labels, on="barcode", how="left", validate="one_to_one")
    data = data.merge(positions, on="barcode", how="left", validate="one_to_one")
    if len(data) != int(record["n_spots"]):
        raise ValueError("Source spot count differs from prior run audit")
    if data[["compartment", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]].isna().any().any():
        raise ValueError("Missing source pathology or coordinates")
    if not np.allclose(data.row, data.pxl_row_in_fullres) or not np.allclose(data.col, data.pxl_col_in_fullres):
        raise ValueError("Source pixel coordinate disagreement")
    return data


def patient_summary(frame: pd.DataFrame, dataset: str, context: str) -> tuple[dict, pd.DataFrame]:
    valid = frame.loc[frame.status.eq("PASS") & frame.dataset.eq(dataset)].copy()
    if context != "all":
        valid = valid.loc[valid.tissue_context.eq(context)]
    patients = valid.groupby("patient_id", as_index=False).pathology_matched_shift_effect_fraction.mean()
    values = patients.pathology_matched_shift_effect_fraction.to_numpy(float)
    if len(values) < 2:
        raise ValueError(f"Only {len(values)} eligible patients in {dataset}/{context}")
    rng = np.random.default_rng(42)
    draws = rng.choice(values, size=(10000, len(values)), replace=True).mean(axis=1)
    lower, upper = np.quantile(draws, [0.025, 0.975])
    leave_one_out = [(values.sum() - value) / (len(values) - 1) for value in values]
    row = {
        "dataset": dataset, "context": context,
        "estimand": "strict tumor/stroma pathology-pair-matched SS-low-source observed-minus-shifted IS contrast",
        "n_eligible_sections": len(valid), "n_eligible_patients": len(values),
        "n_positive_sections": int((valid.pathology_matched_shift_effect_fraction > 0).sum()),
        "n_positive_patients": int((values > 0).sum()),
        "estimate_fraction": float(values.mean()),
        "ci_lower_fraction": float(lower), "ci_upper_fraction": float(upper),
        "leave_one_patient_out_min_fraction": float(min(leave_one_out)),
        "leave_one_patient_out_max_fraction": float(max(leave_one_out)),
        "ci_method": "patient bootstrap percentile, 10000 draws, seed 42",
        "inference_status": "exploratory independent-cohort spatial assessment; translations are structured negative controls",
    }
    patients.insert(0, "dataset", dataset)
    patients.insert(1, "context", context)
    return row, patients


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-zip", type=Path, required=True)
    parser.add_argument("--crosswalk-csv", type=Path, required=True)
    parser.add_argument("--external-qc-json", type=Path, required=True)
    parser.add_argument("--source-run", type=Path, required=True,
                        help="Newly reconstructed 37-section source run audit CSV")
    parser.add_argument("--source-spots", type=Path, required=True,
                        help="Newly reconstructed source spot score tables")
    parser.add_argument("--source-annotations", type=Path, required=True,
                        help="Newly exported author pathology spot labels")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    qc = json.loads(args.external_qc_json.read_text(encoding="utf-8"))
    if len(qc["sections"]) != 12 or len({item["patient_id"] for item in qc["sections"]}) != 6:
        raise ValueError("Unexpected external section/patient denominator")
    if any(item["gate"] != "PASS_METADATA_JOINS_ONLY" for item in qc["sections"]):
        raise ValueError("At least one external input failed effect-blind QC")
    if sha256(args.crosswalk_csv) != qc["crosswalk_sha256"] or sha256(PANEL_SCRIPT) != qc["panel_script_sha256"]:
        raise ValueError("Frozen pathology crosswalk or panel script changed since QC")
    mapping = join_mod.crosswalk(args.crosswalk_csv)
    source_run = pd.read_csv(args.source_run, dtype={"sample": str, "patient_id": str})
    if len(source_run) != 37 or source_run.patient_id.nunique() != 11 or not source_run.status.eq("PASS").all():
        raise ValueError("Unexpected GSE294385 section/patient denominator")
    args.out.mkdir(parents=True, exist_ok=False)
    spot_dir = args.out / "external_spot_scores"
    spot_dir.mkdir()
    scratch = args.out / ".scratch"
    scratch.mkdir()
    offsets = shift_mod.select_offsets()
    sections, shifts = [], []
    with zipfile.ZipFile(args.annotation_zip) as annotations:
        for item in qc["sections"]:
            section_path = Path(item["section_zip"])
            sample = f"{item['patient_id']}_Rep{item['replicate']}"
            print(f"External {sample}", flush=True)
            try:
                if sha256(section_path) != item["section_sha256"]:
                    raise ValueError("Section archive SHA-256 changed after QC")
                data, score_quality = external_data(section_path, annotations, mapping, scratch)
                data.to_csv(spot_dir / f"Table_{sample}_spot_scores.csv")
                result, section_shifts = matched_effect(data, offsets)
                result.update(score_quality)
            except Exception as error:
                result, section_shifts = {"status": "FAIL", "reason": f"{type(error).__name__}: {error}"}, []
            sections.append({
                "dataset": "Valdeolivas_2024", "sample": sample, "patient_id": item["patient_id"],
                "replicate": item["replicate"], "tissue_context": "Primary_CRC", **result
            })
            shifts.extend({"dataset": "Valdeolivas_2024", "sample": sample, "patient_id": item["patient_id"], **record}
                          for record in section_shifts)
            pd.DataFrame(sections).to_csv(args.out / "Table_external_by_section.csv", index=False)
    for _, row in source_run.iterrows():
        gsm = row["sample"]
        print(f"Source {gsm} {row['patient_id']}", flush=True)
        try:
            data = source_data(row, args.source_spots, args.source_annotations)
            result, section_shifts = matched_effect(data, offsets)
        except Exception as error:
            result, section_shifts = {"status": "FAIL", "reason": f"{type(error).__name__}: {error}"}, []
        sections.append({
            "dataset": "GSE294385", "sample": gsm, "patient_id": row["patient_id"],
            "replicate": np.nan, "tissue_context": row["tissue_context"], **result
        })
        shifts.extend({"dataset": "GSE294385", "sample": gsm, "patient_id": row["patient_id"], **record}
                      for record in section_shifts)
        pd.DataFrame(sections).loc[lambda table: table.dataset.eq("GSE294385")].to_csv(
            args.out / "Table_source_harmonized_by_section.csv", index=False
        )
    pd.DataFrame(shifts).to_csv(args.out / "Table_all_attempted_shifts.csv", index=False)
    all_sections = pd.DataFrame(sections)
    all_sections.to_csv(args.out / "Table_all_by_section.csv", index=False)
    summary_rows, patient_rows = [], []
    for dataset, context in [
        ("Valdeolivas_2024", "all"), ("GSE294385", "all"), ("GSE294385", "Primary_CRC")
    ]:
        try:
            summary, patients = patient_summary(all_sections, dataset, context)
            summary_rows.append(summary)
            patient_rows.append(patients)
        except ValueError as error:
            summary_rows.append({"dataset": dataset, "context": context, "status": "NOT_ESTIMABLE", "reason": str(error)})
    pd.DataFrame(summary_rows).to_csv(args.out / "Table_patient_equal_summary.csv", index=False)
    if patient_rows:
        pd.concat(patient_rows, ignore_index=True).to_csv(args.out / "Table_patient_effects.csv", index=False)
    manifest = {
        "protocol": "protocol/CRC_Valdeolivas_external_spatial_gate_2026-10-08.md",
        "final_lock": "protocol/CRC_Valdeolivas_final_analysis_lock_2026-10-08.md",
        "analysis_script_sha256": sha256(Path(__file__)),
        "crosswalk_sha256": sha256(args.crosswalk_csv),
        "input_qc_sha256": sha256(args.external_qc_json),
        "panel_script_sha256": sha256(PANEL_SCRIPT),
        "shift_script_sha256": sha256(SHIFT_SCRIPT),
        "source_run_audit_sha256": sha256(args.source_run),
        "python_numpy_scipy_pandas": {
            "numpy": np.__version__, "scipy": scipy.__version__, "pandas": pd.__version__
        },
        "sampling_unit": "patient; section effects averaged within patient; spots and translations not independent",
        "inference_limits": ["no linked outcome in external cohort", "structured shifts are not exact permutations"],
    }
    (args.out / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(pd.DataFrame(summary_rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
