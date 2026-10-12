# CRC immune–stromal spatial analysis — frozen code candidate 11.10.0-rc7

This package reruns the **current frozen spatial E4 implementation** and its descriptive edge-selection support audit from public input archives, independently recalculates the patient estimand from the newly generated edge ledger, then renders spatial main Figures 2 and 3. It does not read a saved spot-score table, external QC report, section/patient effect, edge-support table, or figure. The analytic unit for the effect summary is the patient: eligible section effects are averaged within patient before the patient-equal mean and 10,000-draw percentile bootstrap (seed 42).

**Status:** local, upload-ready release candidate; it is not yet publicly deposited. No repository URL, DOI, or article acceptance is asserted. The historical E4 executable recorded in the project archive has a different SHA-256 and has not been recovered; a clean rerun of this current implementation does not prove byte-identical execution of that older program. The edge-selection audit is post hoc and descriptive; it does not revise the E4 effect estimator or the original 30/30 edge threshold. The independent arithmetic check starts from the newly rebuilt edge counts; it does not independently reconstruct raw expression or the graph. The package covers this spatial core, S14b/i/p/q and Figures 2–3, not the bulk-survival meta-analysis, ICI audit, single-cell work, or every manuscript figure/supplement.

## Public inputs (not redistributed)

1. GSE294385: the 37 `GSM*_*.tar.gz` Visium sections in `config/GSE294385_37_section_manifest.csv` and the author `visium_meta_after_qc.rds` metadata. The manifest locks section-to-patient and pathology-label denominators and the RDS MD5. Obtain the original archives from [GEO GSE294385](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294385) and the RDS from the [authors' public metadata directory](https://github.com/yliuup/CRC_micromets_ST/tree/main/Meta_data). Do not silently substitute a different RDS revision.
2. Valdeolivas et al.: the twelve selected section ZIPs listed in `config/valdeolivas_12_sections.txt`, plus `Pathology_SpotAnnotations.zip`, from [Zenodo record 7760264](https://zenodo.org/records/7760264). The separate A798015/S7 pair is not analyzed because it lacks the fixed strict tumor/stroma categories. Additional ZIPs in the input directory are ignored.
3. The effect-blind label crosswalk `config/Table_S11g_pathology_crosswalk.csv` is included as an analysis configuration, not as original expression data. The fixed 18-gene IS and 13-gene SS definitions are in `manuscript/CRC_NKT_v11_1/code_release/scripts/05_spatial_neighborhood_analysis.py`. `config/public_input_sha256.json` locks the exact public archives and analysis-configuration bytes used for this candidate; a changed upstream archive fails before any result is written.

Data remain governed by their original repositories' terms; the MIT license in this package applies **only to our code and documentation**. No downloaded count matrix, image, patient-level clinical record, or third-party annotation archive is included.

## Install and run

The tested environment is Windows 11, CPython 3.14.0. Use a new virtual environment and install `requirements-lock-win-py314.txt`; `requirements.txt` lists direct dependencies for a future resolver on other platforms. This is an exact-version lock, not an offline wheel bundle. The machine needs enough disk space for source matrices, plots, and intermediate CSVs.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-lock-win-py314.txt
.\.venv\Scripts\python.exe -m pip check
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe run_from_public.py `
  --gse-archive-dir "D:\path\to\GSE294385" `
  --gse-metadata-rds "D:\path\to\visium_meta_after_qc.rds" `
  --valdeolivas-section-dir "D:\path\to\Valdeolivas_2024" `
  --annotation-zip "D:\path\to\Pathology_SpotAnnotations.zip" `
  --out "D:\path\to\NEW_output_directory"
```

The output directory **must not exist**. The runner fails if required files are missing, their hashes differ from the frozen reference, source RDS MD5 or manifest denominators disagree, any of the 37 source sections fails reconstruction, any of the 12 external sections fails effect-blind barcode/coordinate/QC joins, or the three expected patient summaries are not estimable. It preserves all attempted sections and translations, including failures. The support gate requires 20 candidate sections, 1,920 attempted translations, 15,360 pair-by-IS-group records, 80 S14p section-pair and 44 S14q patient-pair records, 132 S14i section/patient/cohort records, 672/723 accepted source/external offsets, and no qualifying external cross-pathology pair under the fixed 30/30 threshold. The separate `independent_estimand_from_new_ledger/` outputs recompute offset, section, patient, cohort, leave-one-patient-out and bootstrap values from the new pair ledger, with failed pairs treated as missing, not zero. `input_manifest.json` records SHA-256 values; `run_manifest.json` records output hashes, status, and numerical summaries; `figures/figure_manifest.json` records figure inputs and hashes. Ten numbered stage logs identify failures. Generated local logs or intermediate CSVs may contain absolute input paths; review before sharing outputs.

## Frozen analysis choices and interpretation

- IS/SS are within-section means of standardized `log1p(count / spot total UMI × 10,000)` gene expression, dichotomized at each section's 75th percentile. The 18/13 genes are fixed in code. The multiplier is **10,000**, not one million; the older “CPM” shorthand in historical run metadata was inaccurate.
- A six-neighbor spatial graph is distance-trimmed at 1.6 times the section's median nearest-neighbor pitch. Only explicit author pathology categories mapped to tumor or stroma enter the strict analysis.
- The SS-low source contrast uses identical source/target pathology pairs, at least 30 IS-high and 30 IS-low edges per pair, minimum-count weighting, and 96 fixed lattice translations (16 each at hex radius 5–10). A translation needs at least 100 matched IS-high edges, 10% pre-pair coverage, and at least one eligible pair; a section needs at least 24 accepted translations. These offsets are structured negative controls, **not** independent replicates or an exact permutation test.
- The result is an exploratory small spatial association, not direct cell-contact proof, treatment-response prediction, prognostic validation, or a causal effect. Patient inference must not use individual spots or translations as independent observations.
- S14b records all pathology-pair and IS-group edge accounting for attempted offsets. S14p/q summarize pre/post gate support, including failures and NA, without re-estimating effects. S14i gives section/patient/cohort pair estimability; an undefined pair-specific estimate remains NA even when its mixture contribution is zero. The tissue-edge variable is a graph-degree proxy rather than a histologic edge or annotated tumor boundary. Figure 3 is a post hoc selection diagnostic; offsets are technical controls, never patients.

## Software provenance and validation

The `scripts/` programs and fixed score module are the frozen current implementation. `run_from_public.py` is the path-explicit orchestration layer; it does not change the E4 score, graph, translation, 30/30 edge gate, or patient estimator. The original pre-effect gate and final analysis lock are included under `protocol/`; their recorded 2026-10-08 script/output hashes refer to the historical project layout, not to this newly parameterized release script. `SHA256_MANIFEST.json` and the companion ZIP hash identify exact released bytes. A clean-environment verification record accompanies the ZIP; check its status rather than assuming success from source availability alone.

To verify package contents, compare every file in `SHA256_MANIFEST.json` with its local SHA-256. To rerun, use a fresh output directory. An immutable public URL or DOI should be added to the manuscript only after the authors actually deposit this exact ZIP and verify the repository metadata. Cite GEO GSE294385 and Zenodo 7760264 for data provenance as well as the eventual article.
