# CRC spatial support analysis — public release candidate v1.5.4

This directory is the prepared upload set for the public GitHub repository `leekaka1126/crc-nkt-th1-spatial-immune-niche` and its Zenodo archive. It contains the reviewed v1.5.4 code artifact and the S1–S23 derived-data archive, plus a small set of key derived CSVs for browser-readable access.

## Scope

- **Code and author-owned documentation:** MIT License.
- **Author-generated derived CSV/JSON outputs and original result figures:** CC BY 4.0, with attribution to Jiaozhen Hao and Minghui Li.
- **Original GEO/Zenodo matrices:** not redistributed. Obtain them from the original accessions listed in the manuscript and in `DATA_SOURCES_AND_RIGHTS.md` inside the code archive.

The release preserves the locked E4/U1/U2/S19/S20 analyses. S20 is a post hoc unique-spot support audit and does not re-estimate the E4 effect, confidence interval, 30/30 threshold or patient denominator. U2 is a mechanism-specific known-truth stress test, not a general spatial-transcriptomics benchmark or an independent third-party replication.

## Files

- `Additional_file_4_code_release_v1.5.4.zip` — scoped code, protocols, environments, validation scripts and selected ledgers.
- `Additional_file_2_R10_V5_S23_CANDIDATE.zip` — S1–S23 machine-readable derived tables and audit files.
- `derived_tables/` — selected support, U2 and S20 tables for direct browsing.
- `SHA256_MANIFEST.json` — hashes of every file in this prepared set.

## Reproduction boundary

Start with `RUN_AND_VERIFY.md` inside the code archive. The package-only verifiers reproduce the frozen ledgers and denominators. The historical S20 raw-coordinate reconstruction requires the full local workspace and the public raw inputs; the package documents this explicitly. Full per-replicate U2 JSON checkpoints are not redistributed, while the CSV ledgers and scenario summaries are included. No DOI is claimed until the public repository commit and Zenodo record are verified.

## Source datasets

The analysis uses public sources including GEO GSE294385 and GSE283052 and the Valdeolivas Zenodo record (record 7760264). Source-data licensing and attribution remain governed by the original repositories.
