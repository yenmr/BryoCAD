# Dataset manifests

The S1A/S1B/S1E filenames are historical working identifiers, not the final supplementary-table numbering. They describe 126 bryophyte assemblies, 22 benchmark species and 132 angiosperms respectively. Final manuscript tables can replace these after reconciliation.

`training_genomes_17.tsv` is selected from S1A using `training_status == seen species`. Input paths in source manifests refer to the original infrastructure and must be resolved locally. They are not URLs. `site_qc/` contains input-generation counts for seven model organisms; do not substitute them for final per-model scored-site counts in S3.

Known metadata review: Aristolochia fimbriata and Genlisea aurea accession labels were previously inconsistent with manifest filenames. Resolve those two entries against actual assemblies before release. No accession is corrected by guessing.
