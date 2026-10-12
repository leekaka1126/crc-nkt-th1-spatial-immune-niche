#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Section-level NKT-Th1/stroma neighborhood analysis for verified Space Ranger data.

This script treats each Visium section as one inferential unit. Spots are used
to estimate within-section effects and are never pooled as independent patients.
It can always test IS/SS spatial-neighbour enrichment. Tumour/stroma distance
and boundary analyses additionally require explicit, archived spot annotations;
the script will not infer them from the IS/SS scores.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import pandas as pd
from scipy import io, sparse, stats
from scipy.spatial import cKDTree

matplotlib.use("Agg")
import matplotlib.pyplot as plt


IMMUNE_GENES = [
    "NKG7", "KLRC1", "KLRB1", "GZMK", "GZMB", "PRF1", "TBX21", "IFNG",
    "ICOS", "STAT4", "CXCR3", "CXCL13", "IL7R", "IL32", "IFITM1", "CD3E",
    "GIMAP7", "HAVCR2",
]
MATRIX_GENES = [
    "COL1A1", "COL1A2", "COL3A1", "SPARC", "BGN", "LUM", "IGFBP7",
    "FAP", "ACTA2", "DCN", "P4HA1", "LOX", "FN1",
]


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def decode(values: np.ndarray) -> list[str]:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in values]


def read_10x_panel(matrix_path: Path, panel_genes: list[str]) -> pd.DataFrame:
    """Read only panel genes from a 10x HDF5 matrix or Matrix Market directory."""
    if matrix_path.is_file():
        with h5py.File(matrix_path, "r") as handle:
            matrix = handle["matrix"]
            barcodes = decode(matrix["barcodes"][:])
            features = matrix["features"]
            gene_key = "name" if "name" in features else "gene_names"
            genes = pd.Series(decode(features[gene_key][:])).str.upper().tolist()
            counts = sparse.csc_matrix(
                (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
                shape=tuple(matrix["shape"][:]),
            )
    else:
        matrix_file = matrix_path / "matrix.mtx.gz"
        barcode_file = matrix_path / "barcodes.tsv.gz"
        feature_file = matrix_path / "features.tsv.gz"
        required = [matrix_file, barcode_file, feature_file]
        missing_files = [str(path) for path in required if not path.exists()]
        if missing_files:
            raise FileNotFoundError(f"10x Matrix Market directory is incomplete: {missing_files}")
        counts = io.mmread(matrix_file).tocsc()
        barcodes = pd.read_csv(barcode_file, sep="\t", header=None, dtype=str).iloc[:, 0].tolist()
        features = pd.read_csv(feature_file, sep="\t", header=None, dtype=str)
        if features.shape[1] < 2:
            raise RuntimeError("10x features.tsv.gz must have at least gene ID and gene name columns")
        genes = features.iloc[:, 1].str.upper().tolist()
        if counts.shape != (len(genes), len(barcodes)):
            raise RuntimeError("10x matrix shape does not match its feature or barcode files")
    index_by_gene: dict[str, int] = {}
    for idx, gene in enumerate(genes):
        index_by_gene.setdefault(gene, idx)
    missing = sorted(set(panel_genes) - set(index_by_gene))
    if missing:
        raise RuntimeError(f"Matrix lacks panel genes: {missing}")
    selected = [index_by_gene[gene] for gene in panel_genes]
    library_size = np.asarray(counts.sum(axis=0)).ravel()
    library_size[library_size == 0] = 1.0
    dense = counts[selected, :].toarray().T
    log_expr = np.log1p(dense / library_size[:, None] * 1e4)
    return pd.DataFrame(log_expr, index=barcodes, columns=panel_genes)


def read_coordinates(path: Path) -> pd.DataFrame:
    """Support both modern headered and legacy headerless Space Ranger files."""
    if path.suffix.lower() == ".parquet":
        data = pd.read_parquet(path)
    else:
        data = pd.read_csv(path)
        if "barcode" not in data.columns:
            data = pd.read_csv(path, header=None)
            if data.shape[1] != 6:
                raise RuntimeError("Headerless tissue_positions file must have six columns")
            data.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    aliases = {
        "barcode": "barcode",
        "pxl_row_in_fullres": "row",
        "pxl_col_in_fullres": "col",
        "in_tissue": "in_tissue",
    }
    available = {source: target for source, target in aliases.items() if source in data.columns}
    if not {"barcode", "row", "col"}.issubset(set(available.values())):
        raise RuntimeError("Coordinate file lacks barcode, pxl_row_in_fullres, or pxl_col_in_fullres")
    data = data.rename(columns=available)
    if "in_tissue" in data.columns:
        data = data[pd.to_numeric(data["in_tissue"], errors="coerce") == 1]
    data["barcode"] = data["barcode"].astype(str)
    return data.set_index("barcode")[["row", "col"]].apply(pd.to_numeric, errors="coerce").dropna()


def zmean(expression: pd.DataFrame, genes: list[str]) -> pd.Series:
    values = expression[genes]
    scale = values.std(axis=0, ddof=1).replace(0, np.nan)
    return values.sub(values.mean(axis=0), axis=1).div(scale, axis=1).mean(axis=1)


def read_annotations(path: Path, sample: str, barcodes: pd.Index, out: Path) -> pd.DataFrame:
    """Require human/archived compartment labels; write a template if unavailable."""
    if not path.exists():
        template = pd.DataFrame({"barcode": barcodes, "compartment": ""})
        template.to_csv(out / f"Template_{sample}_spot_compartment_annotation.csv", index=False)
        raise RuntimeError("Spot annotation is required. A blank barcode template was written; label tumour/stroma/other and rerun.")
    labels = pd.read_csv(path, dtype=str)
    if "sample" in labels.columns:
        labels = labels[labels["sample"] == sample]
    if not {"barcode", "compartment"}.issubset(labels.columns):
        raise RuntimeError("Annotation file must contain barcode and compartment columns")
    labels["compartment"] = labels["compartment"].str.lower().str.strip()
    allowed = {"tumor", "stroma", "other"}
    invalid = sorted(set(labels["compartment"].dropna()) - allowed)
    if invalid:
        raise RuntimeError(f"Unsupported compartment labels: {invalid}; allowed: {sorted(allowed)}")
    return labels.drop_duplicates("barcode").set_index("barcode")[["compartment"]].reindex(barcodes)


def neighbour_table(coords: np.ndarray, k: int) -> np.ndarray:
    if len(coords) <= k:
        raise RuntimeError(f"Only {len(coords)} tissue spots; k={k} neighbours is invalid")
    return cKDTree(coords).query(coords, k=k + 1)[1][:, 1:]


def permutation_neighbor_enrichment(high_is: np.ndarray, high_ss: np.ndarray, neighbors: np.ndarray, n_perm: int, seed: int) -> tuple[float, float, float]:
    observed = float(high_ss[neighbors][high_is].mean())
    rng = np.random.default_rng(seed)
    permuted = np.empty(n_perm)
    for i in range(n_perm):
        shuffled = rng.permutation(high_ss)
        permuted[i] = shuffled[neighbors][high_is].mean()
    p_value = (1 + np.sum(np.abs(permuted - permuted.mean()) >= abs(observed - permuted.mean()))) / (n_perm + 1)
    return observed, float(permuted.mean()), float(p_value)


def boundary_midpoints(coords: np.ndarray, compartments: np.ndarray, neighbors: np.ndarray) -> np.ndarray:
    points = []
    for i, group in enumerate(compartments):
        if group != "tumor":
            continue
        for j in neighbors[i]:
            if compartments[j] == "stroma":
                points.append((coords[i] + coords[j]) / 2.0)
    return np.asarray(points)


def plot_spatial_scores(
    scores: pd.DataFrame, image_path: Path | None, scalefactors_path: Path | None, output_prefix: Path
) -> None:
    """Create a visual audit of IS/SS spot scores; image overlay is optional."""
    x = scores["col"].to_numpy(float)
    y = scores["row"].to_numpy(float)
    image = None
    if image_path and scalefactors_path:
        if not image_path.exists() or not scalefactors_path.exists():
            raise FileNotFoundError("Both --histology-image and --scalefactors-json must exist for an image overlay")
        image = plt.imread(image_path)
        scales = json.loads(scalefactors_path.read_text(encoding="utf-8"))
        scale_key = "tissue_lowres_scalef" if "lowres" in image_path.name.lower() else "tissue_hires_scalef"
        if scale_key not in scales:
            raise RuntimeError(f"Missing {scale_key} in {scalefactors_path}")
        x = x * float(scales[scale_key])
        y = y * float(scales[scale_key])

    has_compartments = "compartment" in scores.columns and scores["compartment"].notna().any()
    n_panels = 4 if has_compartments else 3
    fig, axes = plt.subplots(1, n_panels, figsize=(5.3 * n_panels, 5.4), constrained_layout=True)
    panels = [("IS_within_section", "IS (within-section z-mean)", "RdBu_r"),
              ("SS_within_section", "SS (within-section z-mean)", "BrBG")]
    for axis, (column, title, cmap) in zip(axes[:2], panels):
        if image is not None:
            axis.imshow(image, origin="upper")
        values = scores[column].to_numpy(float)
        bound = np.nanquantile(np.abs(values), 0.98)
        scatter = axis.scatter(x, y, c=values, s=4, cmap=cmap, vmin=-bound, vmax=bound, linewidths=0, alpha=0.9)
        fig.colorbar(scatter, ax=axis, fraction=0.046, pad=0.03)
        axis.set_title(title)
        axis.set_aspect("equal")
        axis.invert_yaxis() if image is None else None
        axis.set_axis_off()

    categories = np.select(
        [scores["IS_high_q75"] & scores["SS_high_q75"], scores["IS_high_q75"], scores["SS_high_q75"]],
        ["IS-high + SS-high", "IS-high only", "SS-high only"], default="Neither"
    )
    palette = {"Neither": "#d9d9d9", "IS-high only": "#d73027", "SS-high only": "#4575b4", "IS-high + SS-high": "#7b3294"}
    axis = axes[2]
    if image is not None:
        axis.imshow(image, origin="upper")
    for category in ("Neither", "IS-high only", "SS-high only", "IS-high + SS-high"):
        mask = categories == category
        axis.scatter(x[mask], y[mask], s=4, color=palette[category], linewidths=0, alpha=0.85, label=category)
    axis.set_title("Top-quartile spatial pattern")
    axis.set_aspect("equal")
    axis.invert_yaxis() if image is None else None
    axis.set_axis_off()
    axis.legend(loc="lower left", frameon=True, markerscale=2, fontsize=8)
    if has_compartments:
        axis = axes[3]
        if image is not None:
            axis.imshow(image, origin="upper")
        pathology_palette = {"other": "#bdbdbd", "stroma": "#4575b4", "tumor": "#d73027"}
        for compartment in ("other", "stroma", "tumor"):
            mask = scores["compartment"].eq(compartment).to_numpy()
            axis.scatter(x[mask], y[mask], s=4, color=pathology_palette[compartment], linewidths=0, alpha=0.85, label=compartment)
        axis.set_title("Author pathology compartments")
        axis.set_aspect("equal")
        axis.invert_yaxis() if image is None else None
        axis.set_axis_off()
        axis.legend(loc="lower left", frameon=True, markerscale=2, fontsize=8)
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    root = project_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix", type=Path, required=True,
        help="filtered_feature_bc_matrix.h5 or a filtered_feature_bc_matrix directory",
    )
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument(
        "--annotation", type=Path,
        help="Optional CSV with barcode, compartment, and optional sample; required for tumour-stroma boundary analyses",
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--neighbors", type=int, default=6)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--out", type=Path, default=root / "03_Results" / "02_Tables" / "public_enhancement")
    parser.add_argument("--histology-image", type=Path, help="Optional tissue_hires_image.png or tissue_lowres_image.png for figure overlay")
    parser.add_argument("--scalefactors-json", type=Path, help="Required with --histology-image")
    parser.add_argument("--fig-out", type=Path, default=root / "03_Results" / "01_Figures" / "Public_Enhancement")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.fig_out.mkdir(parents=True, exist_ok=True)
    if bool(args.histology_image) != bool(args.scalefactors_json):
        raise RuntimeError("Use --histology-image and --scalefactors-json together")

    panel = IMMUNE_GENES + MATRIX_GENES
    expression = read_10x_panel(args.matrix, panel)
    coords = read_coordinates(args.coordinates)
    common = expression.index.intersection(coords.index)
    if len(common) < 50:
        raise RuntimeError(f"Only {len(common)} matched tissue spots; at least 50 are required")
    expression = expression.loc[common]
    coords = coords.loc[common]
    n_expression_coordinate_spots = len(common)
    n_unannotated_dropped = 0
    labels = None
    if args.annotation:
        labels = read_annotations(args.annotation, args.sample, common, args.out)
        annotated = labels["compartment"].notna()
        n_unannotated_dropped = int((~annotated).sum())
        if n_unannotated_dropped:
            # Do not fabricate a class for an unrepresented barcode. The
            # author-annotated subset is retained as the sole boundary-analysis
            # population, and the exclusion count is reported in the output.
            common = common[annotated.to_numpy()]
            expression = expression.loc[common]
            coords = coords.loc[common]
            labels = labels.loc[common]
        if len(common) < 50:
            raise RuntimeError(f"Only {len(common)} author-annotated tissue spots remain; at least 50 are required")
        if not {"tumor", "stroma"}.issubset(set(labels["compartment"])):
            raise RuntimeError("Both tumour and stroma annotations are required for boundary analysis")
    else:
        pd.DataFrame({"barcode": common, "compartment": ""}).to_csv(
            args.out / f"Template_{args.sample}_spot_compartment_annotation.csv", index=False
        )

    scores = pd.DataFrame(index=common)
    scores["IS_within_section"] = zmean(expression, IMMUNE_GENES)
    scores["SS_within_section"] = zmean(expression, MATRIX_GENES)
    scores = scores.join(coords)
    if labels is not None:
        scores = scores.join(labels)
    else:
        scores["compartment"] = pd.NA
    scores["sample"] = args.sample
    scores["IS_high_q75"] = scores["IS_within_section"] >= scores["IS_within_section"].quantile(0.75)
    scores["SS_high_q75"] = scores["SS_within_section"] >= scores["SS_within_section"].quantile(0.75)

    xy = scores[["row", "col"]].to_numpy(float)
    neighbors = neighbour_table(xy, args.neighbors)
    observed, null_mean, p_value = permutation_neighbor_enrichment(
        scores["IS_high_q75"].to_numpy(bool), scores["SS_high_q75"].to_numpy(bool), neighbors, args.permutations, 20261004
    )
    summary_row = {
        "sample": args.sample,
        "n_spots": len(scores),
        "n_expression_coordinate_spots_before_annotation_filter": n_expression_coordinate_spots,
        "n_unannotated_spots_dropped": n_unannotated_dropped,
        "k_neighbors": args.neighbors,
        "IS_high_neighbor_SS_high_fraction": observed,
        "permutation_null_fraction": null_mean,
        "within_section_permutation_p": p_value,
        "inference_note": "One section-level effect; combine sections by section-level meta-analysis or mixed effects, never spot pooling.",
    }
    if labels is None:
        scores["distance_to_tumor_stroma_boundary_px"] = np.nan
        summary_row.update({
            "n_tumor_spots": np.nan,
            "n_stroma_spots": np.nan,
            "within_section_spot_level_boundary_distance_rank_sum_p_exploratory": np.nan,
            "boundary_distance_statistical_note": "Not estimable without author pathology compartments.",
            "tumor_IS_high_median_boundary_distance_px": np.nan,
            "tumor_IS_low_median_boundary_distance_px": np.nan,
            "analysis_scope": "COORDINATE_NEIGHBORHOOD_ONLY",
            "decision": "Neighbour enrichment completed. Tumour-stroma distance, boundary, and immune-exclusion claims are blocked pending explicit spot compartments; a barcode template was written.",
        })
    else:
        boundaries = boundary_midpoints(xy, scores["compartment"].to_numpy(str), neighbors)
        if len(boundaries) == 0:
            raise RuntimeError("No tumour-stroma adjacency edges at the selected neighborhood size")
        tree = cKDTree(boundaries)
        scores["distance_to_tumor_stroma_boundary_px"] = tree.query(xy)[0]
        tumour = scores[scores["compartment"] == "tumor"].copy()
        distance_test = stats.mannwhitneyu(
            tumour.loc[tumour["IS_high_q75"], "distance_to_tumor_stroma_boundary_px"],
            tumour.loc[~tumour["IS_high_q75"], "distance_to_tumor_stroma_boundary_px"],
            alternative="two-sided",
        )
        summary_row.update({
            "n_tumor_spots": int((scores["compartment"] == "tumor").sum()),
            "n_stroma_spots": int((scores["compartment"] == "stroma").sum()),
            "tumor_IS_high_q75_fraction": scores.loc[scores["compartment"] == "tumor", "IS_high_q75"].mean(),
            "stroma_IS_high_q75_fraction": scores.loc[scores["compartment"] == "stroma", "IS_high_q75"].mean(),
            "tumor_SS_high_q75_fraction": scores.loc[scores["compartment"] == "tumor", "SS_high_q75"].mean(),
            "stroma_SS_high_q75_fraction": scores.loc[scores["compartment"] == "stroma", "SS_high_q75"].mean(),
            "within_section_spot_level_boundary_distance_rank_sum_p_exploratory": distance_test.pvalue,
            "boundary_distance_statistical_note": "Exploratory spot-level rank-sum contrast within one section; spatial autocorrelation means this is not a patient-level or independent-spot p-value.",
            "tumor_IS_high_median_boundary_distance_px": tumour.loc[tumour["IS_high_q75"], "distance_to_tumor_stroma_boundary_px"].median(),
            "tumor_IS_low_median_boundary_distance_px": tumour.loc[~tumour["IS_high_q75"], "distance_to_tumor_stroma_boundary_px"].median(),
            "analysis_scope": "COORDINATE_PLUS_TUMOUR_STROMA_BOUNDARY",
            "decision": "Neighbour and tumour-stroma boundary analyses completed; retain this as one section-level result.",
        })
    summary = pd.DataFrame([summary_row])
    scores.index.name = "barcode"
    scores.to_csv(args.out / f"Table_{args.sample}_spot_scores_and_boundaries.csv")
    summary.to_csv(args.out / f"Table_{args.sample}_section_neighborhood_effect.csv", index=False)
    plot_spatial_scores(
        scores, args.histology_image, args.scalefactors_json,
        args.fig_out / f"Fig_{args.sample}_IS_SS_spatial_maps"
    )
    print(f"Completed one section: {args.sample}; do not interpret this one row as a patient-level validation.")


if __name__ == "__main__":
    main()
