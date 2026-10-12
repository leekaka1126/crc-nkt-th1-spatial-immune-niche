# Valdeolivas CRC spatial cohort: blinded external-validation gate

Created 2026-10-08. Status: **annotation crosswalk locked and the two-section A416371 technical pilot passed; all-other-section QC and the final analysis lock remain pending**. No IS/SS effect from this cohort has been calculated in this workspace. The source-category crosswalk was fixed without reading expression values or effects in `protocol/CRC_Valdeolivas_pathology_crosswalk_lock_2026-10-08.md`; changes after this point must be logged as exploratory.

## Question and independent unit

Can the residual association between the fixed immune-activity (IS) and stromal (SS) programs be reproduced in an independent pathology-annotated CRC Visium cohort? Patients, not spots or technical replicate sections, are the sampling units. The cohort has seven source patient IDs with two sections per ID. S7 is predominantly non-neoplastic in the source paper and may fail a tumor–stroma interface coverage gate; do not count all seven as evaluable tumors in advance.

The existing GSE294385 source cohort contains eight primary-CRC sections from five patients (`manuscript/CRC_NKT_v11_1/supplementary/Table_S2_spatial_section_audit.csv`), nested within the 37-section/11-patient all-site set. Primary CRC and metastasis are not independent strata of different patients and must not be pooled with the new primary-CRC cohort without explicit site-specific reporting.

## Phase 0: inspect without outcome/effect calculation

1. Obtain only `Pathology_SpotAnnotations.zip` (Zenodo v4, expected MD5 `92882eb0db7c2cacc5077cb5fc9578fd`; archive lists 14 per-section CSVs), then one pair of processed Space Ranger section ZIPs. Do not download raw FASTQ/BAM, whole-slide `.ndpi`, intermediary, or deconvolution archives for this gate.
2. Run `scripts/audit_valdeolivas_external_inputs.py` on the annotation ZIP and section ZIPs. Record archive hashes, member names, CSV header and row counts, barcode uniqueness/missingness, observed category values, patient/replicate mapping, and expression/coordinate file presence. The script does not load expression values or compute IS/SS.
3. Build an explicit, manually reviewed category crosswalk to `tumor`, `stroma`, and `other` using author pathology labels alone. Do not use IS/SS, deconvolution, or observed effect direction to choose categories. Exclude necrosis, non-neoplastic mucosa, artifacts, and ambiguous mixtures from the primary tumor–stroma comparison; document every original category.
4. Verify barcode join cardinality and coordinate units on the pilot pair. Check availability of all 18 IS and 13 SS genes in the filtered matrix, tissue-spot coverage, and presence of a tumor–stroma interface. Audit all 14 sections in the same way before declaring the eligible patient denominator. A section without both compartments or sufficient edges is a prespecified QC failure, not a negative result.
5. Freeze the category crosswalk, sample eligibility rules, input file hashes, script commit/hash, estimand, and analysis date **before computing any external score or contrast**. Do not call the current metadata plan a preregistration.

The first archive and label audit is complete: checksum matched; all 14 CSVs, seven patient IDs, and both replicates were present; every one of 39 observed label strings was mapped. On labels alone, 12 sections from six patients contain both strict tumor and strict stroma. Both S7 sections have zero strict tumor and stroma spots. These are **annotation-only candidates**, not final expression/coordinate-qualified sections. The first processed-data pilot will be S3_Col_R / A416371 Rep1 and Rep2 because both have substantial strict tumor and stroma counts. This choice was made without IS/SS scores or spatial effects. Processed ZIPs total 288.8 MB for the pilot pair and 1,274.0 MB for all 12 annotation-only candidate sections, based on the Zenodo manifest's rounded file sizes.

The A416371 pilot is now complete at the **metadata/join/geometry** level (`reports/CRC_Valdeolivas_A416371_pilot_archive_audit_2026-10-08.json`; `reports/CRC_Valdeolivas_A416371_pilot_spot_join_QC_2026-10-08.json`). Both processed ZIP MD5 values match Zenodo. Rep1 has 2,317 filtered spots and Rep2 has 1,803; in each section the filtered, pathology and in-tissue coordinate barcodes agree one-to-one with zero orphans. The H5 and Matrix Market dimensions and gene/barcode order agree; all 18 IS and 13 SS panel symbols occur in both filtered feature lists. Array parity, position uniqueness and hires-image bounds pass. A six-nearest-neighbor graph built on **all** filtered spots before applying strict pathology labels contains 240 Rep1 and 177 Rep2 directed tumor-to-stroma edges after the fixed 1.6-pitch cutoff. These are feasibility counts, not patient-independent observations or evidence of an IS/SS association. Expression counts, detection rates, scores, null translations and spatial effects were not read or computed. The next ten candidate processed ZIPs (five patients; 985.2 MB from rounded manifest sizes) are listed in `reports/CRC_Valdeolivas_remaining_candidate_downloads_2026-10-08.csv`. S7's two ZIPs are not needed for the strict primary comparison.

## Proposed fixed analysis after Phase 0 lock

- Use the existing, unchanged 18-gene IS and 13-gene SS panels in `manuscript/CRC_NKT_v11_1/code_release/scripts/05_spatial_neighborhood_analysis.py`.
- Per spot, use `log1p(count / total_UMI * 10000)`, within-section gene z-scores (`ddof=1`), panel mean, and within-section top-quartile IS-high/SS-high flags. Report gene detection and score distributions; quartiles fix each section's high fraction and cannot establish absolute between-patient expression differences.
- Recreate the six-nearest-neighbor graph; discard edges longer than 1.6 times the section's median nearest-spot pitch. Restrict source spots to SS-low. Match contrasts within harmonized source–target pathology pairs, using only tumor/stroma compartments. Compare observed IS-high-minus-IS-low SS-high neighbor fractions against the same eligible edges under the existing 96 fixed intact-map translations (seed 42). These translations are structured negative controls, **not** an exact permutation test.
- Preserve the existing minimums: at least 30 matched edges in each IS group, 100 matched IS-high edges, 10% matched-edge coverage, and 24 accepted translations. If a minimum is impossible because broad categories differ in granularity from GSE294385 Layer3, stop, record the failure, and amend the protocol *before* any effect computation; do not silently relax thresholds.
- Recompute GSE294385 using the **same harmonized tumor/stroma pair definition**. Its earlier Layer3-matched +1.65 percentage-point result is a different estimand and must not be used as the numeric external comparator.
- Average valid sections within each patient, then weight patients equally within each cohort. Report every section/patient, eligible-denominator flow, patient-equal effect, 10,000 patient-bootstrap percentile CI (seed 42), positive-patient count, and leave-one-patient-out range. With approximately six to seven external patients, emphasize effect size and uncertainty rather than a binary P-value.
- Primary analysis: all eligible harmonized tumor/stroma pathology pairs. Prespecified secondary analyses: directed tumor-to-stroma boundary edges, primary CRC subset in GSE294385, and a physical-distance immune-penetration description if image/coordinate scaling is verified. The secondary analyses do not substitute for a failed primary result.

## Interpretation gates

The result supports external spatial reproducibility only if the annotation/coordinate joins pass, the same estimand is calculable in both cohorts, and patient-level direction and effect size are not dominated by one individual. A nonmatching or negative result is reportable and rules out a straightforward validation claim. No spatial patient here has linked ICI response or recurrence endpoints; do not call this predictive or prognostic validation. No iNKT or Th1 cell identity is inferred from Visium spots.

## Source receipts

- Zenodo v4 data manifest and per-file MD5: https://zenodo.org/records/7760264
- Source paper: https://doi.org/10.1038/s41698-023-00488-4
- Prior pathology-free partial-replication lock: `protocol/GSE283052_partial_spatial_replication_locked_2026-10-08.md`

