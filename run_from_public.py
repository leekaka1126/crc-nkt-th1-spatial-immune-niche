#!/usr/bin/env python
"""Rebuild spatial E4, S14 support, independent estimand audit, Figures 2-3.

No saved spot-score, QC, effect, support, or figure table is an input.
The output directory must not exist. A nonzero exit means the release gate failed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parent
PANEL = ROOT / "manuscript/CRC_NKT_v11_1/code_release/scripts/05_spatial_neighborhood_analysis.py"
MANIFEST = ROOT / "config/GSE294385_37_section_manifest.csv"
CROSSWALK = ROOT / "config/Table_S11g_pathology_crosswalk.csv"
EXTERNAL_LIST = ROOT / "config/valdeolivas_12_sections.txt"
INPUT_REFERENCE = ROOT / "config/public_input_sha256.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def execute(arguments: list[str], log: Path) -> None:
    with log.open("w", encoding="utf-8") as output:
        completed = subprocess.run(arguments, cwd=ROOT, stdout=output,
                                   stderr=subprocess.STDOUT, check=False)
    if completed.returncode:
        raise RuntimeError(f"Stage failed (exit {completed.returncode}); inspect {log}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gse-archive-dir", type=Path, required=True)
    parser.add_argument("--gse-metadata-rds", type=Path, required=True)
    parser.add_argument("--valdeolivas-section-dir", type=Path, required=True)
    parser.add_argument("--annotation-zip", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("--out must be a new directory; archived results are never overwritten")
    samples = pd.read_csv(MANIFEST, dtype=str)
    if len(samples) != 37 or samples.geo_sample_accession.nunique() != 37:
        raise ValueError("Frozen GSE294385 manifest is not 37 unique sections")
    gse_archives = [args.gse_archive_dir / f"{row.geo_sample_accession}_{row.author_sample_id}.tar.gz"
                    for row in samples.itertuples(index=False)]
    names = [line.strip() for line in EXTERNAL_LIST.read_text(encoding="utf-8").splitlines()
             if line.strip() and not line.startswith("#")]
    if len(names) != 12 or len(set(names)) != 12:
        raise ValueError("Frozen external section list is not 12 unique ZIP names")
    sections = [args.valdeolivas_section_dir / name for name in names]
    if len(sections) != 12:
        raise ValueError(f"Expected exactly 12 selected Valdeolivas section ZIPs, got {len(sections)}")
    required = [*gse_archives, args.gse_metadata_rds, *sections, args.annotation_zip,
                MANIFEST, CROSSWALK, EXTERNAL_LIST, PANEL]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing inputs: {missing}")
    observed_hashes = {path.name: sha256(path) for path in required}
    expected_hashes = json.loads(INPUT_REFERENCE.read_text(encoding="utf-8"))["input_sha256"]
    missing_reference = sorted(set(expected_hashes) - set(observed_hashes))
    changed = sorted(name for name in expected_hashes.keys() & observed_hashes.keys()
                     if expected_hashes[name] != observed_hashes[name])
    if missing_reference or changed:
        raise ValueError(f"Frozen public input/configuration hash mismatch: "
                         f"missing={missing_reference}, changed={changed}")
    args.out.mkdir(parents=True, exist_ok=False)
    input_manifest = {
        "release": "11.10.0-rc7", "python": sys.version,
        "input_sha256": observed_hashes,
        "frozen_input_reference_sha256": sha256(INPUT_REFERENCE),
        "input_paths_are_local_only": True,
        "source_archive_count": 37, "external_section_count": 12,
    }
    (args.out / "input_manifest.json").write_text(
        json.dumps(input_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    source = args.out / "source_from_raw"
    execute([sys.executable, str(ROOT / "scripts/gse294385_full_spatial.py"),
             "--manifest", str(MANIFEST), "--archive-dir", str(args.gse_archive_dir),
             "--metadata", str(args.gse_metadata_rds), "--analysis-script", str(PANEL),
             "--out", str(source)], args.out / "01_source.log")
    source_table = pd.read_csv(source / "Table_GSE294385_full_spatial_run_audit.csv")
    if len(source_table) != 37 or not source_table.status.eq("PASS").all():
        raise AssertionError("Source reconstruction did not pass all 37 sections; see 01_source.log")

    qc = args.out / "external_input_qc.json"
    command = [sys.executable, str(ROOT / "scripts/audit_valdeolivas_spot_joins.py"),
               "--annotation-zip", str(args.annotation_zip), "--crosswalk-csv", str(CROSSWALK),
               "--panel-script", str(PANEL), "--output-json", str(qc)]
    for section in sections:
        command.extend(("--section-zip", str(section)))
    execute(command, args.out / "02_external_qc.log")
    external_qc = json.loads(qc.read_text(encoding="utf-8"))
    if len(external_qc["sections"]) != 12 or any(
            item["gate"] != "PASS_METADATA_JOINS_ONLY" for item in external_qc["sections"]):
        raise AssertionError("External effect-blind QC did not pass all 12 sections")

    spatial = args.out / "e4_from_raw"
    execute([sys.executable, str(ROOT / "scripts/run_valdeolivas_pathology_matched_replication.py"),
             "--annotation-zip", str(args.annotation_zip), "--crosswalk-csv", str(CROSSWALK),
             "--external-qc-json", str(qc),
             "--source-run", str(source / "Table_GSE294385_full_spatial_run_audit.csv"),
             "--source-spots", str(source / "sections"),
             "--source-annotations", str(source / "annotations"),
             "--out", str(spatial)], args.out / "03_spatial_e4.log")
    effects = pd.read_csv(spatial / "Table_patient_equal_summary.csv")
    expected = {("Valdeolivas_2024", "all"), ("GSE294385", "all"),
                ("GSE294385", "Primary_CRC")}
    if set(zip(effects.dataset, effects.context)) != expected or effects.estimate_fraction.isna().any():
        raise AssertionError("E4 patient summaries are missing or not estimable")

    support = args.out / "edge_support_from_raw"
    execute([sys.executable, str(ROOT / "scripts/audit_v11_8_full_edge_selection.py"),
             "--e4-dir", str(spatial), "--source-dir", str(source),
             "--out", str(support)], args.out / "04_edge_audit.log")
    execute([sys.executable, str(ROOT / "scripts/summarize_v11_8_full_edge_selection.py"),
             "--audit-dir", str(support), "--e4-dir", str(spatial)],
            args.out / "05_edge_summary.log")
    execute([sys.executable, str(ROOT / "scripts/build_crc_v11_10_R2k_S14_support_rollup.py"),
             "--audit-dir", str(support)], args.out / "06_support_rollup.log")
    audit_manifest = json.loads((support / "manifest.json").read_text(encoding="utf-8"))
    rollup_manifest = json.loads((support / "Table_S14p-q_support_rollup_manifest.json").read_text(encoding="utf-8"))
    if (audit_manifest["candidate_sections"], audit_manifest["attempted_shifts"],
            audit_manifest["pair_by_IS_group_rows"], audit_manifest["frozen_attempts_verified"]) != (20, 1920, 15360, 1920):
        raise AssertionError("Full edge audit failed its frozen denominator gate")
    rollup = rollup_manifest["summary"]
    if (rollup["section_rows"], rollup["patient_rows"]) != (80, 44):
        raise AssertionError("S14p/q section or patient denominators differ")
    pair_checks = {(item["dataset"], item["pathology_pair"]): item
                   for item in rollup["cohort_pair_checks"]}
    for dataset, accepted in (("GSE294385", 672), ("Valdeolivas_2024", 723)):
        for pair in ("tumor->tumor", "tumor->stroma", "stroma->tumor", "stroma->stroma"):
            if pair_checks[(dataset, pair)]["analysis_offset_n"] != accepted:
                raise AssertionError(f"Accepted-offset denominator differs for {dataset} {pair}")
    if any(pair_checks[("Valdeolivas_2024", pair)]["pair_qualified_analysis_offset_n"]
           for pair in ("tumor->stroma", "stroma->tumor")):
        raise AssertionError("External cross-pathology support is not zero")

    execute([sys.executable, str(ROOT / "scripts/build_v11_10_revision1_pair_ledger.py"),
             "--audit-dir", str(support)], args.out / "07_multilevel_estimability.log")
    estimability = pd.read_csv(support / "Table_S14i_multilevel_pair_estimability.csv")
    if (len(estimability), estimability.scope.eq("section").sum(),
            estimability.scope.eq("patient").sum(),
            estimability.scope.eq("cohort").sum()) != (132, 80, 44, 8):
        raise AssertionError("S14i multilevel row denominator differs")
    external_cross = estimability.loc[
        estimability.scope.eq("cohort") & estimability.dataset.eq("Valdeolivas_2024") &
        estimability.pathology_pair.isin(("tumor->stroma", "stroma->tumor"))]
    if (len(external_cross) != 2 or external_cross.pair_eligible.any() or
            external_cross.estimate.notna().any() or
            not external_cross.mixture_weight.eq(0).all()):
        raise AssertionError("S14i external cross-pathology NA/zero distinction failed")

    independent = args.out / "independent_estimand_from_new_ledger"
    execute([sys.executable, str(ROOT / "scripts/check_crc_R2k4_patient_estimand_independent.py"),
             "--pair-ledger", str(support / "Table_shift_pair_and_edge_hops.csv"),
             "--e4-shifts", str(spatial / "Table_all_attempted_shifts.csv"),
             "--e4-sections", str(spatial / "Table_all_by_section.csv"),
             "--e4-patients", str(spatial / "Table_patient_effects.csv"),
             "--e4-summary", str(spatial / "Table_patient_equal_summary.csv"),
             "--out", str(independent)], args.out / "08_independent_estimand.log")
    independent_report = json.loads((independent / "audit_report.json").read_text(encoding="utf-8"))
    if (independent_report["status"], independent_report["n_pair_group_rows"],
            independent_report["n_offsets"]) != ("PASS", 15360, 1920):
        raise AssertionError("Independent patient-estimand arithmetic failed")
    independent_cohorts = {(item["dataset"], item["context"]): item
                           for item in independent_report["cohorts"]}
    if (independent_cohorts[("Valdeolivas_2024", "all")]["n_accepted_section_offsets"],
            independent_cohorts[("GSE294385", "Primary_CRC")]["n_accepted_section_offsets"]) != (723, 672):
        raise AssertionError("Independent estimator accepted-offset denominator differs")

    figures = args.out / "figures"
    execute([sys.executable, str(ROOT / "scripts/plot_valdeolivas_patient_effects.py"),
             "--run-dir", str(spatial), "--output-prefix", str(figures / "Fig2_patient_contrasts")],
            args.out / "09_figure2.log")
    execute([sys.executable, str(ROOT / "scripts/make_v11_10_edge_selection_figure.py"),
             "--audit-dir", str(support), "--figure-dir", str(figures)],
            args.out / "10_figure3.log")
    figure_paths = [figures / f"{stem}.{ext}" for stem in
                    ("Fig2_patient_contrasts", "Fig3_v11_10_edge_selection_audit")
                    for ext in ("pdf", "png")]
    if any(not path.is_file() or path.stat().st_size == 0 for path in figure_paths):
        raise AssertionError("Figure 2 or Figure 3 PDF/PNG was not produced")
    figure_manifest = {
        "role": "spatial_main_figures_only",
        "figures": [
            {"figure": "Figure 2", "stem": "Fig2_patient_contrasts",
             "source_tables": ["e4_from_raw/Table_patient_effects.csv",
                               "e4_from_raw/Table_patient_equal_summary.csv"]},
            {"figure": "Figure 3", "stem": "Fig3_v11_10_edge_selection_audit",
             "source_tables": ["edge_support_from_raw/Table_patient_equal_cohort_selection.csv",
                               "edge_support_from_raw/Table_shift_edge_selection.csv",
                               "edge_support_from_raw/Table_shift_pair_and_edge_hops.csv"]}],
        "output_sha256": {path.name: sha256(path) for path in figure_paths},
        "formats": "PDF vector and PNG 610 dpi",
    }
    (figures / "figure_manifest.json").write_text(
        json.dumps(figure_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (figures / "_figure_manifest.md").write_text(
        "# Spatial figure manifest\n\n"
        "Both figures are generated from this run, not copied from manuscript assets. "
        "Journal-specific dimensions remain to be checked at submission.\n\n"
        "| Figure | Files | Analysis source | Role |\n"
        "| --- | --- | --- | --- |\n"
        "| 2 | `Fig2_patient_contrasts.pdf`, `.png` | "
        "`e4_from_raw/Table_patient_effects.csv`; "
        "`e4_from_raw/Table_patient_equal_summary.csv` | Patient-level E4 contrast |\n"
        "| 3 | `Fig3_v11_10_edge_selection_audit.pdf`, `.png` | "
        "`edge_support_from_raw/` audit tables | Post hoc support diagnostic |\n\n"
        "Exact SHA-256 values are in `figure_manifest.json`. PNG export is 610 dpi; "
        "PDF output is vector graphics.\n", encoding="utf-8")
    output_manifest = {
        "status": "PASS_RAW_TO_E4_S14BIPQ_INDEPENDENT_ESTIMAND_FIG2_FIG3", "release": "11.10.0-rc7",
        "n_source_sections": 37, "n_external_sections": 12,
        "patient_summary": effects.to_dict(orient="records"),
        "edge_audit": {key: audit_manifest[key] for key in
                       ("candidate_sections", "attempted_shifts", "pair_by_IS_group_rows",
                        "frozen_attempts_verified")},
        "support_rollup": rollup,
        "s14i_rows": {"total": len(estimability), "section": 80, "patient": 44, "cohort": 8},
        "independent_estimand": independent_report["cohorts"],
        "output_sha256": {path.relative_to(args.out).as_posix(): sha256(path)
                          for path in sorted(args.out.rglob("*.csv"))},
        "figure_sha256": figure_manifest["output_sha256"],
        "historical_e4_executable_recovered": False,
        "warning": "Current frozen code rerun; not byte-identical proof of the unavailable historical E4 executable",
    }
    (args.out / "run_manifest.json").write_text(
        json.dumps(output_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(effects.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
