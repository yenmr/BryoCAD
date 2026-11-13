# Changelog

## Unreleased

- Organize continued-pretraining, motif, variant, accessibility and gene-annotation workflows.
- Provide FASTA/GFF3 representative-transcript protein reconstruction.
- Include the GeneCAD v0.4.0 modification patch and published model links.
- Document the 123-assembly selected-contig evaluation and UTR-independent coding-transcript metrics.
- Add examples and selected figure-source tables; further toy datasets are in preparation.

- Simplified ATAC training to a direct 4,000-step run.
- Standardized public script names and updated references and checksums.

- Fixed SNP VCF list exhaustion and short-contig centering; added regression tests and score-header declaration.
- Added Apache-2.0 license, upstream license copies and manuscript-based citation metadata.
- Added four real-sequence motif examples and the released ATAC LoRA model link.
- Simplified protein-preparation documentation to the FASTA/GFF3 workflow.

- Added published annotation downloads, verified public HF revisions and ATAC task-head presence.
- Renamed figure-specific analysis scripts and added Linux environment recipes and executable CPU/GPU release checks.

- Use FP32 for V100 smoke tests; preserve explicit precision overrides and capture native crash diagnostics.
- Split zero-shot, ATAC and GeneCAD runners, with automatic HF downloads and an aggregate failure summary.

- Restore custom-forward compatibility in the ATAC smoke test; verify saved LoRA and classification-head tensors after loading.

- Added raw TSV input and metadata-preserving TSV output to the original LoRA predict CLI; added explicit synthetic example labels and TSV checks.

- Correct replica dispatch in LoRA forward wrappers; restrict ATAC smoke tests to the first allocated GPU and record visible device count.

- Select GeneCAD’s own environment independently; preflight Python requirements, imports and CUDA before downloads; limit the smoke-test pipeline to one attempt.

- Record environment packages without requiring pip; verify an actual basic CUDA operation instead of relying only on device visibility.

- Add a unified conda build based on the successful dicos environment, including native architecture patches, pinned runtime packages, a dicos build job and environment capture. GPU validation of the clean unified build remains pending.
