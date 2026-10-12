#!/usr/bin/env python
"""Effect-blind expression/annotation/coordinate QC for Valdeolivas Visium ZIPs.

Reads matrix metadata and panel gene names, but never reads count values or
computes IS, SS, spatial contrasts, or outcomes. Input ZIPs are not modified.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import struct
import sys
import zipfile

import h5py
import numpy as np
from scipy.spatial import cKDTree


SAMPLE_RE = re.compile(r"_(A\d+)_Rep([12])(?:_X)?\.zip$", re.IGNORECASE)
BARCODE_RE = re.compile(r"^[ACGT]{16}-\d+$")


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def exactly_one(archive: zipfile.ZipFile, suffix: str) -> str:
    matches = [name for name in archive.namelist() if name.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"Expected one {suffix!r} in {archive.filename}; found {matches}")
    return matches[0]


def gz_lines(archive: zipfile.ZipFile, member: str):
    with archive.open(member) as raw, gzip.GzipFile(fileobj=raw) as decoded:
        for line in io.TextIOWrapper(decoded, encoding="utf-8", newline=""):
            yield line.rstrip("\r\n")


def panel_from_script(path: Path) -> dict[str, list[str]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    requested = {"IMMUNE_GENES", "MATRIX_GENES"}
    found = {}
    for statement in tree.body:
        if isinstance(statement, ast.Assign):
            for target in statement.targets:
                if isinstance(target, ast.Name) and target.id in requested:
                    found[target.id] = ast.literal_eval(statement.value)
    if set(found) != requested:
        raise ValueError(f"Panel assignments absent in {path}")
    return found


def annotation_for_sample(archive: zipfile.ZipFile, patient_id: str, rep: int) -> tuple[str, dict[str, str]]:
    pattern = re.compile(rf"_{patient_id}_Rep{rep}(?:_X)?\.csv$", re.IGNORECASE)
    members = [name for name in archive.namelist() if pattern.search(name)]
    if len(members) != 1:
        raise ValueError(f"Expected one annotation CSV for {patient_id} Rep{rep}: {members}")
    member = members[0]
    with archive.open(member) as raw, io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        barcode_fields = [field for field in fields if field.strip().lower() == "barcode"]
        if len(barcode_fields) != 1 or len(fields) != 2:
            raise ValueError(f"Unexpected annotation columns in {member}: {fields}")
        label_field = next(field for field in fields if field != barcode_fields[0])
        pairs = [(row[barcode_fields[0]].strip(), row[label_field].strip()) for row in reader]
    barcodes = [barcode for barcode, _ in pairs]
    if len(set(barcodes)) != len(barcodes) or any(not BARCODE_RE.fullmatch(value) for value in barcodes):
        raise ValueError(f"Nonunique or malformed annotation barcode in {member}")
    return member, dict(pairs)


def crosswalk(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    labels = [row["original_label"] for row in rows]
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate pathology label in crosswalk")
    mapping = {row["original_label"]: row["primary_class"] for row in rows}
    if set(mapping.values()) - {"tumor", "stroma", "other"}:
        raise ValueError("Unexpected crosswalk class")
    return mapping


def matrix_metadata(archive: zipfile.ZipFile) -> dict:
    barcodes_member = exactly_one(archive, "/filtered_feature_bc_matrix/barcodes.tsv.gz")
    features_member = exactly_one(archive, "/filtered_feature_bc_matrix/features.tsv.gz")
    matrix_member = exactly_one(archive, "/filtered_feature_bc_matrix/matrix.mtx.gz")
    h5_member = exactly_one(archive, "/filtered_feature_bc_matrix.h5")
    tsv_barcodes = [line.strip() for line in gz_lines(archive, barcodes_member)]
    tsv_features = [line.split("\t") for line in gz_lines(archive, features_member)]
    if any(len(parts) < 2 for parts in tsv_features):
        raise ValueError("Filtered features TSV has fewer than two columns")
    with archive.open(matrix_member) as raw, gzip.GzipFile(fileobj=raw) as decoded:
        text = io.TextIOWrapper(decoded, encoding="ascii")
        header = text.readline().strip()
        dimensions = next(line.strip() for line in text if line.strip() and not line.startswith("%"))
    if not header.startswith("%%MatrixMarket matrix coordinate"):
        raise ValueError(f"Unexpected Matrix Market header: {header}")
    mm_shape = [int(item) for item in dimensions.split()]
    if len(mm_shape) != 3:
        raise ValueError(f"Unexpected Matrix Market dimensions: {dimensions}")
    with h5py.File(io.BytesIO(archive.read(h5_member)), "r") as handle:
        matrix = handle["matrix"]
        h5_shape = [int(value) for value in matrix["shape"][:]]
        h5_barcodes = [value.decode() if isinstance(value, bytes) else str(value) for value in matrix["barcodes"][:]]
        names = matrix["features"]["name"][:]
        h5_genes = [value.decode() if isinstance(value, bytes) else str(value) for value in names]
        stored_nonzero = int(matrix["data"].shape[0])
    report = {
        "filtered_barcodes": len(tsv_barcodes),
        "filtered_features": len(tsv_features),
        "matrix_market_shape_and_nnz": mm_shape,
        "h5_shape": h5_shape,
        "h5_stored_nonzero": stored_nonzero,
        "tsv_h5_barcodes_identical_order": tsv_barcodes == h5_barcodes,
        "tsv_h5_gene_names_identical_order": [parts[1] for parts in tsv_features] == h5_genes,
        "barcode_duplicates": len(tsv_barcodes) - len(set(tsv_barcodes)),
        "barcode_format_failures": sum(not BARCODE_RE.fullmatch(value) for value in tsv_barcodes),
    }
    report["shape_consistent"] = (
        h5_shape == mm_shape[:2] == [len(tsv_features), len(tsv_barcodes)]
        and stored_nonzero == mm_shape[2]
    )
    return report, tsv_barcodes, [parts[1] for parts in tsv_features]


def position_metadata(archive: zipfile.ZipFile) -> tuple[dict, dict[str, tuple[int, int, int, float, float]]]:
    member = exactly_one(archive, "/spatial/tissue_positions_list.csv")
    with archive.open(member) as raw, io.TextIOWrapper(raw, encoding="utf-8", newline="") as text:
        rows = list(csv.reader(text))
    if any(len(row) != 6 for row in rows):
        raise ValueError(f"Unexpected coordinate row width in {member}")
    positions = {}
    duplicates = 0
    for barcode, in_tissue, array_row, array_col, pixel_row, pixel_col in rows:
        if barcode in positions:
            duplicates += 1
        positions[barcode] = (int(in_tissue), int(array_row), int(array_col), float(pixel_row), float(pixel_col))
    tissue = {barcode: spot for barcode, spot in positions.items() if spot[0] == 1}
    lattice = [(spot[1], spot[2]) for spot in positions.values()]
    scales_member = exactly_one(archive, "/spatial/scalefactors_json.json")
    scales = json.loads(archive.read(scales_member))
    image_member = exactly_one(archive, "/spatial/tissue_hires_image.png")
    with archive.open(image_member) as raw:
        png_header = raw.read(24)
    if png_header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"Invalid PNG header: {image_member}")
    image_width, image_height = struct.unpack(">II", png_header[16:24])
    hires_scale = float(scales["tissue_hires_scalef"])
    image_oob = sum(
        not (0 <= spot[3] * hires_scale < image_height and 0 <= spot[4] * hires_scale < image_width)
        for spot in tissue.values()
    )
    return {
        "all_coordinate_rows": len(rows),
        "in_tissue_rows": len(tissue),
        "coordinate_duplicate_barcodes": duplicates,
        "coordinate_duplicate_array_cells": len(lattice) - len(set(lattice)),
        "invalid_tissue_flags": sum(spot[0] not in {0, 1} for spot in positions.values()),
        "invalid_visium_lattice_parity": sum((spot[1] + spot[2]) % 2 != 0 for spot in positions.values()),
        "hires_scale": hires_scale,
        "hires_image_width_height": [image_width, image_height],
        "in_tissue_centers_outside_hires_image": image_oob,
    }, positions


def graph_screen(barcodes: list[str], positions: dict, labels: dict, mapping: dict) -> dict:
    xy = np.asarray([(positions[barcode][3], positions[barcode][4]) for barcode in barcodes], dtype=float)
    distance, neighbor = cKDTree(xy).query(xy, k=7)
    pitch = float(np.median(distance[:, 1]))
    accepted = distance[:, 1:] <= 1.6 * pitch
    source = np.repeat(np.arange(len(barcodes)), 6)[accepted.reshape(-1)]
    target = neighbor[:, 1:].reshape(-1)[accepted.reshape(-1)]
    categories = [mapping[labels[barcode]] for barcode in barcodes]
    pair_counts = Counter(f"{categories[i]}->{categories[j]}" for i, j in zip(source, target))
    strict_counts = {key: pair_counts.get(key, 0) for key in ("tumor->tumor", "tumor->stroma", "stroma->tumor", "stroma->stroma")}
    return {
        "graph_spots": len(barcodes),
        "median_nearest_neighbor_pitch_fullres_px": round(pitch, 6),
        "max_edge_distance_fullres_px": round(1.6 * pitch, 6),
        "accepted_directed_edges": len(source),
        "strict_directed_edge_counts": strict_counts,
        "strict_directed_edges_total": sum(strict_counts.values()),
        "method": "6-nearest neighbors on all filtered spots; retain edges <=1.6 median nearest-neighbor pitch; classify only after graph construction",
    }


def one_section(section_path: Path, annotation_zip: zipfile.ZipFile, mapping: dict, panels: dict) -> dict:
    match = SAMPLE_RE.search(section_path.name)
    if not match:
        raise ValueError(f"Cannot parse patient and replicate from {section_path.name}")
    patient_id, rep = match.group(1), int(match.group(2))
    annotation_member, labels = annotation_for_sample(annotation_zip, patient_id, rep)
    with zipfile.ZipFile(section_path) as archive:
        matrix_report, barcodes, genes = matrix_metadata(archive)
        position_report, positions = position_metadata(archive)
    matrix_set, annotation_set, position_set = set(barcodes), set(labels), set(positions)
    tissue_set = {barcode for barcode, spot in positions.items() if spot[0] == 1}
    missing_labels = sorted(set(labels.values()) - set(mapping))
    class_counts = Counter(mapping.get(labels[barcode], "unmapped") for barcode in matrix_set & annotation_set)
    gene_set = {gene.upper() for gene in genes}
    panel_report = {
        panel: {"required": len(required), "present": len(set(required) & gene_set), "missing": sorted(set(required) - gene_set)}
        for panel, required in panels.items()
    }
    join = {
        "annotation_member": annotation_member,
        "annotation_barcodes": len(annotation_set),
        "matrix_barcodes": len(matrix_set),
        "coordinate_barcodes": len(position_set),
        "tissue_coordinate_barcodes": len(tissue_set),
        "matrix_without_annotation": len(matrix_set - annotation_set),
        "annotation_without_matrix": len(annotation_set - matrix_set),
        "matrix_without_coordinates": len(matrix_set - position_set),
        "matrix_flagged_outside_tissue": len(matrix_set - tissue_set),
        "joined_class_counts": dict(class_counts),
        "unmapped_annotation_labels": missing_labels,
    }
    report = {
        "section_zip": str(section_path.resolve()),
        "section_sha256": sha256(section_path),
        "patient_id": patient_id,
        "replicate": rep,
        "matrix": matrix_report,
        "positions": position_report,
        "join": join,
        "panels": panel_report,
        "graph": None,
        "issues": [],
        "not_checked": ["expression count values", "gene detection rates", "IS/SS scores", "spatial effect estimates"],
    }
    if not matrix_report["shape_consistent"] or not matrix_report["tsv_h5_barcodes_identical_order"] or not matrix_report["tsv_h5_gene_names_identical_order"]:
        report["issues"].append("Filtered TSV, Matrix Market, and H5 metadata disagree")
    if matrix_report["barcode_duplicates"] or matrix_report["barcode_format_failures"]:
        report["issues"].append("Filtered matrix barcodes are duplicated or malformed")
    if any(position_report[key] for key in (
        "coordinate_duplicate_barcodes", "coordinate_duplicate_array_cells", "invalid_tissue_flags",
        "invalid_visium_lattice_parity", "in_tissue_centers_outside_hires_image"
    )):
        report["issues"].append("Coordinate, lattice, or image geometry integrity issue")
    if any(join[key] for key in (
        "matrix_without_annotation", "annotation_without_matrix", "matrix_without_coordinates",
        "matrix_flagged_outside_tissue", "unmapped_annotation_labels"
    )):
        report["issues"].append("Expression, annotation, and coordinate barcodes do not form a complete one-to-one tissue join")
    if any(panel["missing"] for panel in panel_report.values()):
        report["issues"].append("At least one frozen IS/SS gene symbol is absent")
    if not report["issues"]:
        report["graph"] = graph_screen(barcodes, positions, labels, mapping)
    report["gate"] = "PASS_METADATA_JOINS_ONLY" if not report["issues"] else "HOLD_FOR_REVIEW"
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-zip", type=Path, required=True)
    parser.add_argument("--crosswalk-csv", type=Path, required=True)
    parser.add_argument("--panel-script", type=Path, required=True)
    parser.add_argument("--section-zip", type=Path, action="append", required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    for path in [args.annotation_zip, args.crosswalk_csv, args.panel_script, *args.section_zip]:
        if not path.is_file():
            parser.error(f"Not a file: {path}")
    try:
        mapping = crosswalk(args.crosswalk_csv)
        panels = panel_from_script(args.panel_script)
        with zipfile.ZipFile(args.annotation_zip) as annotation_zip:
            sections = [one_section(path, annotation_zip, mapping, panels) for path in args.section_zip]
        result = {
            "purpose": "effect-blind pilot input and geometric QC",
            "crosswalk_sha256": sha256(args.crosswalk_csv),
            "panel_script_sha256": sha256(args.panel_script),
            "sections": sections,
        }
    except (OSError, ValueError, KeyError, UnicodeError, zipfile.BadZipFile, csv.Error, RuntimeError) as error:
        print(f"QC failed: {error}", file=sys.stderr)
        return 2
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0 if all(item["gate"] == "PASS_METADATA_JOINS_ONLY" for item in sections) else 1


if __name__ == "__main__":
    raise SystemExit(main())
