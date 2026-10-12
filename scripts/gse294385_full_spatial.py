#!/usr/bin/env python
"""Reproduce a patient-aware GSE294385 spatial validation from public archives.

The score, q75 cut and section-level algorithm are imported unchanged from the
existing project script. Only explicit author pathology labels are used.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import shutil
import subprocess
import sys
import tarfile
import tempfile
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
from scipy.spatial import cKDTree


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def export_author_annotations(metadata: pd.DataFrame, author_id: str, gsm: str, dest: Path) -> tuple[Path, dict]:
    prefix = author_id + "~"
    subset = metadata.loc[metadata["Barcode"].str.startswith(prefix), ["Barcode", "Layer1", "Layer2", "Layer3"]].copy()
    if subset.empty:
        raise ValueError(f"No author metadata rows for {author_id}")
    subset["barcode"] = subset["Barcode"].str[len(prefix):]
    if subset["barcode"].duplicated().any():
        raise ValueError(f"Duplicate author barcodes for {author_id}")
    layer = subset["Layer3"].fillna("")
    is_tumor = layer.str.contains(r"^(?:Liver|Lung) (?:macro|micro)metastasis tumor$|^Tumor cells ", case=False, regex=True)
    is_stroma = layer.str.contains(r"^(?:Liver|Lung) (?:macro|micro)metastasis stroma$|^Tumor stroma ", case=False, regex=True)
    subset["compartment"] = np.select([is_tumor, is_stroma], ["tumor", "stroma"], default="other")
    subset["sample"] = gsm
    result = subset[["barcode", "compartment", "sample", "Layer1", "Layer2", "Layer3"]].sort_values("barcode")
    path = dest / f"Table_{gsm}_author_pathology_spot_compartments.csv"
    result.to_csv(path, index=False)
    counts = result["compartment"].value_counts().to_dict()
    return path, {"n_author_spots": len(result), "n_author_tumor": counts.get("tumor", 0), "n_author_stroma": counts.get("stroma", 0)}


def extract_required(archive_path: Path, author_id: str, temp_dir: Path) -> tuple[Path, Path, dict]:
    with tarfile.open(archive_path, "r:gz") as archive:
        members = {m.name.replace("\\", "/"): m for m in archive.getmembers() if m.isfile()}
        prefix = author_id + "/"
        if not members or any(not name.startswith(prefix) for name in members):
            raise ValueError("Archive has unexpected top-level path or no files")
        h5_name = prefix + "filtered_feature_bc_matrix.h5"
        mtx_names = [prefix + "filtered_feature_bc_matrix/" + name for name in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")]
        if h5_name in members:
            needed = [h5_name]
            matrix_path = temp_dir / "filtered_feature_bc_matrix.h5"
        elif all(name in members for name in mtx_names):
            needed = list(mtx_names)
            matrix_path = temp_dir / "filtered_feature_bc_matrix"
            matrix_path.mkdir()
        else:
            raise ValueError("No complete filtered 10x expression matrix")
        coord_candidates = [prefix + "spatial/" + name for name in ("tissue_positions.csv", "tissue_positions_list.csv", "tissue_positions.parquet")]
        coord_name = next((name for name in coord_candidates if name in members), None)
        if not coord_name:
            raise ValueError("No spot coordinate file")
        needed.append(coord_name)
        for name in needed:
            member = members[name]
            if member.size <= 0:
                raise ValueError(f"Empty required member: {name}")
            if name == h5_name:
                target = matrix_path
            elif name in mtx_names:
                target = matrix_path / Path(name).name
            else:
                target = temp_dir / Path(name).name
            with archive.extractfile(member) as source, target.open("wb") as output:
                shutil.copyfileobj(source, output, length=2**20)
            if target.stat().st_size != member.size:
                raise IOError(f"Extracted size mismatch: {name}")
    return matrix_path, temp_dir / Path(coord_name).name, {"n_archive_files": len(members), "matrix_member": h5_name if h5_name in members else " | ".join(mtx_names), "coordinate_member": coord_name}


def patient_summary(sections: pd.DataFrame, output: Path) -> None:
    good = sections.loc[sections["status"] == "PASS"].copy()
    if good.empty:
        return
    good["neighbor_excess_fraction"] = good["IS_high_neighbor_SS_high_fraction"] - good["permutation_null_fraction"]
    good["stroma_minus_tumor_IS_high_fraction"] = good["stroma_IS_high_q75_fraction"] - good["tumor_IS_high_q75_fraction"]
    good.to_csv(output / "Table_GSE294385_37_section_effects.csv", index=False)
    patient = good.groupby("patient_id", as_index=False).agg(
        n_sections=("sample", "size"),
        neighbor_excess_fraction=("neighbor_excess_fraction", "mean"),
        stroma_minus_tumor_IS_high_fraction=("stroma_minus_tumor_IS_high_fraction", "mean"),
        boundary_distance_contrast_spot_units=("boundary_distance_contrast_spot_units", "mean"),
    )
    patient.to_csv(output / "Table_GSE294385_11_patient_spatial_effects.csv", index=False)
    rng = np.random.default_rng(42)
    rows = []
    for label, column in [
        ("IS-high to SS-high six-neighbor excess, fraction", "neighbor_excess_fraction"),
        ("Stroma minus tumor IS-high prevalence, fraction", "stroma_minus_tumor_IS_high_fraction"),
        ("Tumor IS-high minus IS-low boundary distance, spot units", "boundary_distance_contrast_spot_units"),
    ]:
        values = patient[column].dropna().to_numpy(float)
        n = len(values)
        if not n:
            continue
        resamples = rng.choice(values, size=(10000, n), replace=True).mean(axis=1)
        ci = np.quantile(resamples, [0.025, 0.975])
        if n <= 20:
            signs = np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int8)
            null = (signs * values).mean(axis=1)
            p = float((np.abs(null) >= abs(values.mean()) - 1e-12).mean())
        else:
            p = np.nan
        rows.append({"estimand": label, "unit": "patient", "n_patients": n,
                     "estimate": float(values.mean()), "ci_lower": float(ci[0]), "ci_upper": float(ci[1]),
                     "sign_flip_two_sided_p_exploratory": p,
                     "ci_method": "patient bootstrap percentile, 10000 resamples, seed 42",
                     "test_note": "Exploratory patient-level sign-flip; section permutation does not preserve spatial autocorrelation"})
    pd.DataFrame(rows).to_csv(output / "Table_GSE294385_patient_level_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--analysis-script", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0, help="Pilot count; 0 runs all 37")
    parser.add_argument("--permutations", type=int, default=999)
    args = parser.parse_args()
    for folder in (args.out, args.out / "annotations", args.out / "sections", args.out / "figures", args.out / ".scratch"):
        folder.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, dtype={"geo_sample_accession": str, "author_sample_id": str, "patient_id": str})
    required = {"geo_sample_accession", "author_sample_id", "patient_id", "archive_bytes", "n_author_annotated_spots", "n_author_tumor_spots", "n_author_stroma_spots"}
    if not required.issubset(manifest.columns) or manifest["patient_id"].isna().any():
        raise ValueError("Manifest lacks required columns or patient mapping")
    if len(manifest) != 37 or manifest["geo_sample_accession"].duplicated().any():
        raise ValueError("Manifest must contain 37 unique GSMs")
    metadata = pyreadr.read_r(str(args.metadata))[None]
    if not {"Barcode", "Layer1", "Layer2", "Layer3"}.issubset(metadata.columns):
        raise ValueError("Author RDS has unexpected columns")
    source_md5 = md5(args.metadata)
    expected_md5 = str(manifest["source_metadata_md5"].iloc[0])
    if source_md5 != expected_md5:
        raise ValueError(f"RDS MD5 mismatch: {source_md5} != {expected_md5}")
    run_rows = []
    for i, row in enumerate(manifest.itertuples(index=False), 1):
        if args.limit and i > args.limit:
            break
        gsm, author = row.geo_sample_accession, row.author_sample_id
        archive_path = args.archive_dir / f"{gsm}_{author}.tar.gz"
        rec = {"sample": gsm, "author_sample_id": author, "patient_id": row.patient_id,
               "tissue_context": row.tissue_context, "archive": str(archive_path), "status": "FAIL"}
        try:
            if not archive_path.is_file() or archive_path.stat().st_size != int(row.archive_bytes):
                raise ValueError("Archive absent or byte size differs from GEO manifest")
            annotation, counts = export_author_annotations(metadata, author, gsm, args.out / "annotations")
            rec.update(counts)
            for observed, expected, name in [
                (counts["n_author_spots"], row.n_author_annotated_spots, "author spots"),
                (counts["n_author_tumor"], row.n_author_tumor_spots, "tumor spots"),
                (counts["n_author_stroma"], row.n_author_stroma_spots, "stroma spots"),
            ]:
                if observed != int(expected):
                    raise ValueError(f"{name} count differs from locked manifest: {observed} != {expected}")
            with tempfile.TemporaryDirectory(prefix=f"{gsm}_", dir=args.out / ".scratch") as temp:
                matrix, coords, archive_info = extract_required(archive_path, author, Path(temp))
                rec.update(archive_info)
                if counts["n_author_tumor"] == 0 or counts["n_author_stroma"] == 0:
                    raise ValueError("No explicit tumor and stroma spots in author annotation")
                command = [sys.executable, str(args.analysis_script), "--matrix", str(matrix),
                           "--coordinates", str(coords), "--annotation", str(annotation),
                           "--sample", gsm, "--permutations", str(args.permutations),
                           "--out", str(args.out / "sections"), "--fig-out", str(args.out / "figures")]
                completed = subprocess.run(command, capture_output=True, text=True)
                if completed.returncode:
                    raise RuntimeError((completed.stderr or completed.stdout)[-2500:])
            section = pd.read_csv(args.out / "sections" / f"Table_{gsm}_section_neighborhood_effect.csv").iloc[0].to_dict()
            rec.update(section)
            spots = pd.read_csv(args.out / "sections" / f"Table_{gsm}_spot_scores_and_boundaries.csv")
            xy = spots[["row", "col"]].to_numpy(float)
            nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
            spacing = float(np.median(nn[nn > 0]))
            if spacing <= 0 or not np.isfinite(spacing):
                raise ValueError("No finite nearest-neighbor spot spacing")
            rec["median_nearest_spot_spacing_px"] = spacing
            rec["boundary_distance_contrast_spot_units"] = (
                rec["tumor_IS_high_median_boundary_distance_px"] -
                rec["tumor_IS_low_median_boundary_distance_px"]
            ) / spacing
            rec["status"] = "PASS"
            rec["error"] = ""
        except Exception as error:
            rec["error"] = f"{type(error).__name__}: {error}"[-3000:]
            print(traceback.format_exc(limit=2), file=sys.stderr, flush=True)
        run_rows.append(rec)
        pd.DataFrame(run_rows).to_csv(args.out / "Table_GSE294385_full_spatial_run_audit.csv", index=False)
        print(f"[{i}/{min(len(manifest), args.limit or len(manifest))}] {gsm} {rec['status']}: {rec.get('error', '')[:160]}", flush=True)
    sections = pd.DataFrame(run_rows)
    patient_summary(sections, args.out)
    (args.out / "run_manifest.json").write_text(json.dumps({
        "metadata_rds_md5": source_md5,
        "analysis_script": str(args.analysis_script), "analysis_script_md5": md5(args.analysis_script),
        "runner_md5": md5(Path(__file__)), "python": sys.version,
        "permutations": args.permutations, "neighbors": 6,
        "score": "Existing 18-gene IS and 13-gene SS log1p(count/spot_total_UMI*10000) within-section z-mean; q75 per section",
        "interpretation": "37 sections nested in 11 patients; metastatic CRC spatial evidence, not prognostic or ICI outcome validation",
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"FINAL: {int((sections.status == 'PASS').sum())} passed / {len(sections)} attempted; {sections.loc[sections.status == 'PASS', 'patient_id'].nunique()} patients", flush=True)


if __name__ == "__main__":
    main()
