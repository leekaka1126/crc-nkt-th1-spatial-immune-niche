#!/usr/bin/env python
"""Post hoc, exact-reproduction audit of pathology-pair and tissue-edge selection.

No score or effect rule is changed. Shifts and directed edges are not inference units.
"""
from __future__ import annotations

import argparse
from collections import deque
import hashlib
import json
from pathlib import Path
import platform
import warnings

import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SCRIPT = ROOT / "scripts/run_valdeolivas_pathology_matched_replication.py"
PAIR_NAMES = ("tumor->tumor", "tumor->stroma", "stroma->tumor", "stroma->stroma")
BIN_NAMES = ("0", "1", "2", "3plus")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tissue_edge_hops(n_spots: int, all_source: np.ndarray, all_target: np.ndarray) -> np.ndarray:
    """Shortest 6-NN graph hops to a sampled-tissue-footprint boundary proxy."""
    adjacency = [set() for _ in range(n_spots)]
    for source, target in zip(all_source, all_target):
        source, target = int(source), int(target)
        if source != target:
            adjacency[source].add(target)
            adjacency[target].add(source)
    boundary = [spot for spot, neighbors in enumerate(adjacency) if len(neighbors) < 6]
    if not boundary:
        raise ValueError("No degree-<6 tissue-footprint boundary spots")
    hops = np.full(n_spots, -1, dtype=np.int32)
    hops[boundary] = 0
    queue = deque(boundary)
    while queue:
        source = queue.popleft()
        for target in adjacency[source]:
            if hops[target] < 0:
                hops[target] = hops[source] + 1
                queue.append(target)
    if (hops < 0).any():
        raise ValueError("Graph component lacks sampled-tissue-footprint boundary")
    return hops


def counts_by_pair_group(mask: np.ndarray, code: np.ndarray) -> np.ndarray:
    return np.bincount(code[mask], minlength=8).astype(np.int64)


def four_bins(mask: np.ndarray, hops: np.ndarray) -> tuple[int, int, int, int]:
    if not mask.any():
        return (0, 0, 0, 0)
    values = np.minimum(hops[mask], 3)
    result = np.bincount(values, minlength=4)
    return tuple(int(value) for value in result)


def store_bins(row: dict, label: str, mask: np.ndarray, hops: np.ndarray) -> None:
    for name, value in zip(BIN_NAMES, four_bins(mask, hops)):
        row[f"{label}_hop_{name}"] = value


def one_section(data: pd.DataFrame, identity: dict, offsets: list[tuple[int, int, int]],
                ref: pd.DataFrame, frozen_section: pd.Series) -> tuple[list[dict], list[dict], dict]:
    original = __import__("importlib.util").util
    spec = original.spec_from_file_location("frozen_pathology_matched", SOURCE_SCRIPT)
    assert spec is not None and spec.loader is not None
    module = original.module_from_spec(spec)
    spec.loader.exec_module(module)
    layer = data.compartment.astype(str).to_numpy()
    immune = module.true_flags(data.IS_high_q75, "IS_high_q75")
    stroma = module.true_flags(data.SS_high_q75, "SS_high_q75")
    ar, ac = data.array_row.to_numpy(int), data.array_col.to_numpy(int)
    if ((ar + ac) % 2).any() or data[["array_row", "array_col"]].duplicated().any():
        raise ValueError(f"Bad Visium lattice: {identity}")
    distance, neighbor = cKDTree(data[["row", "col"]].to_numpy(float)).query(
        data[["row", "col"]].to_numpy(float), k=7)
    pitch = float(np.median(distance[:, 1]))
    valid = (distance[:, 1:] <= 1.6 * pitch).reshape(-1)
    all_source = np.repeat(np.arange(len(data)), 6)[valid]
    all_target = neighbor[:, 1:].reshape(-1)[valid]
    strict = np.isin(layer[all_source], ["tumor", "stroma"]) & np.isin(
        layer[all_target], ["tumor", "stroma"])
    source, target = all_source[strict], all_target[strict]
    if len(source) < 100 or int(immune[source].sum()) < 100:
        raise ValueError(f"Section fails pre-shift edge gate: {identity}")
    if min(int((~stroma[source] & immune[source]).sum()),
           int((~stroma[source] & ~immune[source]).sum())) < 30:
        raise ValueError(f"Section fails pre-shift SS-low edge gate: {identity}")
    hops = tissue_edge_hops(len(data), all_source, all_target)
    source_hops = hops[source]
    target_hops = hops[target]
    base_eligible = ~stroma[source]
    source_high = immune[source]
    n_selected = int(source_high.sum())
    n_base_high = int((base_eligible & source_high).sum())
    pair = (layer[source] == "stroma").astype(np.int8) * 2 + (layer[target] == "stroma").astype(np.int8)
    code = pair * 2 + source_high.astype(np.int8)
    original_counts = counts_by_pair_group(base_eligible, code)
    baseline, _, _, _ = module.shift_mod.matched_layer_contrast(source, target, immune, stroma, layer)
    lookup = np.full((int(ar.max()) + 1, int(ac.max()) + 1), -1, dtype=np.int32)
    lookup[ar, ac] = np.arange(len(data), dtype=np.int32)
    if int((lookup >= 0).sum()) != len(data):
        raise ValueError("Lattice lookup loses spots")
    shift_rows, pair_rows = [], []
    for dr, dc, radius in offsets:
        rr, cc = ar[target] + dr, ac[target] + dc
        in_grid = (rr >= 0) & (rr < lookup.shape[0]) & (cc >= 0) & (cc < lookup.shape[1])
        mapped = np.full(len(target), -1, dtype=np.int32)
        mapped[in_grid] = lookup[rr[in_grid], cc[in_grid]]
        exists = mapped >= 0
        compatible = exists.copy()
        compatible[exists] &= layer[mapped[exists]] == layer[target[exists]]
        mapped_hops = np.full(len(target), -1, dtype=np.int32)
        mapped_hops[exists] = hops[mapped[exists]]
        n_matched = int((compatible & source_high).sum())
        before_coverage = n_matched / n_selected
        pre_gate = n_matched >= 100 and before_coverage >= 0.10
        compatible_eligible = compatible & base_eligible
        grouped = counts_by_pair_group(compatible_eligible, code)
        pair_pass = (grouped[0::2] >= 30) & (grouped[1::2] >= 30)
        eligible_count = int(compatible_eligible.sum())
        if not pre_gate or eligible_count < 100:
            pair_pass[:] = False
        retained = compatible_eligible & pair_pass[pair]
        final_counts = counts_by_pair_group(retained, code)
        n_pairs = int(pair_pass.sum())
        final_high = int(final_counts[1::2].sum())
        after_coverage = final_high / n_base_high if pre_gate else 0.0
        if not pre_gate:
            status = "FAIL_COVERAGE"
        elif not n_pairs:
            status = "FAIL_PAIR_MINIMUM"
        else:
            status = "PASS"
        pair_deltas, pair_observed, pair_shifted, weights = [], [], [], []
        mapped_eligible = exists & base_eligible
        mapped_counts = counts_by_pair_group(mapped_eligible, code)
        for pair_number, pair_name in enumerate(PAIR_NAMES):
            observations = []
            for high in (0, 1):
                cell = pair_number * 2 + high
                cell_original = base_eligible & (code == cell)
                cell_mapped = mapped_eligible & (code == cell)
                cell_compatible = compatible_eligible & (code == cell)
                cell_retained = retained & (code == cell)
                row = {**identity, "shift_array_row": dr, "shift_array_col": dc,
                       "hex_radius": radius, "status": status, "pathology_pair": pair_name,
                       "IS_high_source": bool(high), "n_original_SS_low_edges": int(original_counts[cell]),
                       "n_mapped_in_tissue": int(mapped_counts[cell]),
                       "n_off_tissue_or_lattice": int(original_counts[cell] - mapped_counts[cell]),
                       "n_pathology_mismatch": int(mapped_counts[cell] - grouped[cell]),
                       "n_compatible": int(grouped[cell]), "pair_passed_minimum": bool(pair_pass[pair_number]),
                       "n_pair_retained": int(final_counts[cell]),
                       "n_observed_SS_high_on_retained": int(stroma[target[cell_retained]].sum()),
                       "n_shifted_SS_high_on_retained": int(stroma[mapped[cell_retained]].sum())}
                store_bins(row, "original_target_original", cell_original, target_hops)
                store_bins(row, "source_original", cell_original, source_hops)
                store_bins(row, "mapped_target_mapped", cell_mapped, mapped_hops)
                store_bins(row, "compatible_original_target", cell_compatible, target_hops)
                store_bins(row, "compatible_source", cell_compatible, source_hops)
                store_bins(row, "compatible_mapped_target", cell_compatible, mapped_hops)
                store_bins(row, "retained_original_target", cell_retained, target_hops)
                store_bins(row, "retained_source", cell_retained, source_hops)
                store_bins(row, "retained_mapped_target", cell_retained, mapped_hops)
                pair_rows.append(row)
                observations.append(row)
            low, high = observations
            if pair_pass[pair_number]:
                observed = high["n_observed_SS_high_on_retained"] / high["n_pair_retained"] - (
                    low["n_observed_SS_high_on_retained"] / low["n_pair_retained"])
                shifted = high["n_shifted_SS_high_on_retained"] / high["n_pair_retained"] - (
                    low["n_shifted_SS_high_on_retained"] / low["n_pair_retained"])
                weight = min(high["n_pair_retained"], low["n_pair_retained"])
                pair_observed.append(observed)
                pair_shifted.append(shifted)
                pair_deltas.append(observed - shifted)
                weights.append(weight)
        observed_contrast = float(np.average(pair_observed, weights=weights)) if weights else np.nan
        shifted_contrast = float(np.average(pair_shifted, weights=weights)) if weights else np.nan
        effect = float(np.average(pair_deltas, weights=weights)) if weights else np.nan
        shift_row = {**identity, "shift_array_row": dr, "shift_array_col": dc, "hex_radius": radius,
                     "status": status, "n_spots": len(data), "n_tissue_boundary_proxy_spots": int((hops == 0).sum()),
                     "n_all_short_directed_edges": len(all_source), "n_strict_short_directed_edges": len(source),
                     "n_original_IS_high_strict_edges": n_selected,
                     "n_original_SS_low_source_edges": int(base_eligible.sum()),
                     "n_original_SS_low_IS_high_edges": n_base_high,
                     "n_mapped_in_tissue_SS_low_source_edges": int(mapped_eligible.sum()),
                     "n_compatible_SS_low_source_edges": eligible_count,
                     "n_pair_retained_SS_low_source_edges": int(retained.sum()),
                     "n_matched_IS_high_edges_before_pair_filter": n_matched,
                     "matched_IS_high_edge_fraction_before_pair_filter": before_coverage,
                     "n_matched_pathology_pairs": n_pairs,
                     "matched_IS_high_edge_fraction_after_pair_filter": after_coverage,
                     "n_pair_retained_IS_high_edges": final_high,
                     "baseline_full_edge_observed_IS_contrast": baseline,
                     "observed_matched_edge_IS_contrast": observed_contrast,
                     "shifted_matched_edge_IS_contrast": shifted_contrast,
                     "observed_matched_minus_full_edge_contrast": observed_contrast - baseline,
                     "observed_minus_shifted_IS_contrast": effect}
        for high in (0, 1):
            suffix = "high" if high else "low"
            mask = base_eligible & (source_high == high)
            store_bins(shift_row, f"original_{suffix}_original_target", mask, target_hops)
            store_bins(shift_row, f"original_{suffix}_source", mask, source_hops)
            store_bins(shift_row, f"retained_{suffix}_original_target", mask & retained, target_hops)
            store_bins(shift_row, f"retained_{suffix}_source", mask & retained, source_hops)
            store_bins(shift_row, f"retained_{suffix}_mapped_target", mask & retained, mapped_hops)
        shift_rows.append(shift_row)
        check = ref.loc[(ref.shift_array_row == dr) & (ref.shift_array_col == dc)]
        if len(check) != 1:
            raise AssertionError(f"Missing or repeated frozen offset in {identity}: {(dr, dc)}")
        check = check.iloc[0]
        for field in ("status", "n_matched_IS_high_edges_before_pair_filter", "n_matched_pathology_pairs"):
            if shift_row[field] != check[field]:
                raise AssertionError(f"Frozen mismatch {identity} {(dr, dc)} {field}: {shift_row[field]} vs {check[field]}")
        for field in ("matched_IS_high_edge_fraction_before_pair_filter",
                      "matched_IS_high_edge_fraction_after_pair_filter", "observed_minus_shifted_IS_contrast"):
            if not np.isclose(shift_row[field], check[field], rtol=1e-10, atol=1e-12, equal_nan=True):
                raise AssertionError(f"Frozen mismatch {identity} {(dr, dc)} {field}: {shift_row[field]} vs {check[field]}")
    accepted = sum(row["status"] == "PASS" for row in shift_rows)
    if len(shift_rows) != 96 or len(pair_rows) != 96 * 8 or accepted != int(frozen_section.n_accepted_shifts):
        raise AssertionError(f"Frozen section count mismatch: {identity}")
    section_row = {**identity, "section_status": frozen_section.status, "n_spots": len(data),
                   "n_tissue_boundary_proxy_spots": int((hops == 0).sum()),
                   "n_strict_short_directed_edges": len(source),
                   "n_original_SS_low_source_edges": int(base_eligible.sum()),
                   "n_accepted_shifts": accepted, "baseline_full_edge_observed_IS_contrast": baseline}
    return shift_rows, pair_rows, section_row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--e4-dir", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    reference_sections = args.e4_dir / "Table_all_by_section.csv"
    reference_shifts = args.e4_dir / "Table_all_attempted_shifts.csv"
    external_dir = args.e4_dir / "external_spot_scores"
    source_run_path = args.source_dir / "Table_GSE294385_full_spatial_run_audit.csv"
    import importlib.util
    spec = importlib.util.spec_from_file_location("frozen_pathology_matched_loader", SOURCE_SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    offsets = module.shift_mod.select_offsets()
    sections = pd.read_csv(reference_sections)
    sections = sections.loc[sections.dataset.eq("Valdeolivas_2024") | (
        sections.dataset.eq("GSE294385") & sections.tissue_context.eq("Primary_CRC"))].copy()
    frozen = pd.read_csv(reference_shifts)
    frozen = frozen.loc[frozen.dataset.eq("Valdeolivas_2024") | (
        frozen.dataset.eq("GSE294385") & frozen["sample"].isin(
            sections.loc[sections.dataset.eq("GSE294385"), "sample"]))].copy()
    if len(sections) != 20 or len(frozen) != 1920 or len(offsets) != 96:
        raise AssertionError("Expected 20 primary candidate sections and 1,920 frozen shifts")
    source_run = pd.read_csv(source_run_path, dtype={"sample": str, "patient_id": str})
    all_shift, all_pair, all_section = [], [], []
    inputs = {path.name: sha256(path) for path in
              [reference_sections, reference_shifts, SOURCE_SCRIPT,
               ROOT / "scripts/gse294385_spatial_autocorr_audit.py", source_run_path]}
    for section in sections.itertuples(index=False):
        identity = {"dataset": section.dataset, "sample": section.sample, "patient_id": section.patient_id}
        section_ref = frozen.loc[frozen.dataset.eq(section.dataset) & frozen["sample"].eq(section.sample)]
        if len(section_ref) != 96:
            raise AssertionError(f"Wrong frozen shift count: {identity}")
        if section.dataset == "GSE294385":
            source_record = source_run.loc[source_run["sample"].eq(section.sample)]
            if len(source_record) != 1:
                raise AssertionError(f"Missing source section: {identity}")
            source_record = source_record.iloc[0]
            data = module.source_data(source_record, args.source_dir / "sections",
                                      args.source_dir / "annotations")
            archive_file = Path(source_record["archive"])
            inputs[str(archive_file)] = sha256(archive_file)
            score_file = args.source_dir / "sections" / f"Table_{section.sample}_spot_scores_and_boundaries.csv"
            label_file = args.source_dir / "annotations" / f"Table_{section.sample}_author_pathology_spot_compartments.csv"
            inputs[score_file.name] = sha256(score_file)
            inputs[label_file.name] = sha256(label_file)
        else:
            score_file = external_dir / f"Table_{section.sample}_spot_scores.csv"
            inputs[score_file.name] = sha256(score_file)
            data = pd.read_csv(score_file)
        print(f"Auditing {section.dataset} {section.sample}", flush=True)
        shifts, pairs, section_row = one_section(data, identity, offsets, section_ref, pd.Series(section._asdict()))
        all_shift.extend(shifts)
        all_pair.extend(pairs)
        all_section.append(section_row)
    args.out.mkdir(parents=True, exist_ok=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        pd.DataFrame(all_shift).to_csv(args.out / "Table_shift_edge_selection.csv", index=False)
        pd.DataFrame(all_pair).to_csv(args.out / "Table_shift_pair_and_edge_hops.csv", index=False)
        pd.DataFrame(all_section).to_csv(args.out / "Table_section_baselines.csv", index=False)
    manifest = {"analysis_role": "post_hoc_descriptive", "candidate_sections": len(all_section),
                "attempted_shifts": len(all_shift), "pair_by_IS_group_rows": len(all_pair),
                "frozen_attempts_verified": len(all_shift), "tissue_edge_proxy": "graph hops to degree <6 spot",
                "no_random_resampling": True, "python": platform.python_version(), "numpy": np.__version__,
                "pandas": pd.__version__, "scipy": scipy.__version__, "inputs_sha256": inputs,
                "script_sha256": sha256(Path(__file__))}
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Verified {len(all_shift)} shifts, {len(all_pair)} pathology-pair/IS-group cells", flush=True)


if __name__ == "__main__":
    main()
