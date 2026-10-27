# Upload and archive instructions

1. Upload this directory to the existing public GitHub repository under `release/v1.5.4/`.
2. Commit the complete directory without renaming the two ZIP archives.
3. Create a Git tag named `v1.5.4` on that commit.
4. Enable the repository in Zenodo and archive the tag. Record the resulting DOI only after an unauthenticated download and SHA-256 check succeed.
5. Update the manuscript Data Availability and Code Availability statements with the verified GitHub commit/tag and Zenodo DOI. Do not insert a placeholder DOI.

The two authors should verify `SHA256_MANIFEST.json` before the commit. Do not upload raw matrices, patient-identifying files, access tokens, local path logs, or the Word manuscript in this release directory.
