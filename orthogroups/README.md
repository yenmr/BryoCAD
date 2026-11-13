# Proteins from genome FASTA and GFF3

Select one coding transcript per gene from the GFF3 annotation, then reconstruct and translate its CDS using the matching genome FASTA.

```bash
python orthogroups/extract_proteins_from_gff.py \
  --fasta MpTak_v6.fa --gff MpTak_v6.1r2.gff.gz \
  --output results/proteomes/M_polymorpha.fa \
  --report results/protein_qc/M_polymorpha.selection.tsv
```

Python standard library only; `.gz` inputs are supported. The FASTA must match the annotation assembly and contig identifiers. Run separately for each assembly; OrthoFinder input files should contain only protein FASTA files, not selection reports.

## Selection and translation

- Default `--length-metric transcript`: longest merged exonic length among CDS-containing transcripts. If exons are absent, use the CDS+UTR union; absent UTR annotations cannot reconstruct missing transcript length.
- Ties: greater CDS length, then lexicographic transcript ID. `--length-metric cds` explicitly selects by CDS length instead; it is not the default.
- Selection happens before translation/QC. An invalid selected transcript is reported; the program does not silently substitute a shorter isoform.
- CDSs are assembled in transcript orientation. Only initial CDS phase is trimmed; subsequent phases are checked and their bases retained to complete split codons. Standard nuclear genetic code (table 1); alternative genetic codes/translation exceptions are not handled.
- A terminal stop is removed; internal stops, invalid phases and out-of-bounds CDSs are rejected with reasons. Partial terminal codons are omitted and reported. Ambiguous codons become X.
- Existing protein QC defaults are retained: minimum 50 aa, X fraction <=0.30, initial X run <20. Adjust only if the actual intended proteome uses different filters.
- FASTA identifiers are GFF gene IDs. Use `--id-suffix '|Mpo'` only when needed to match original orthogroup identifiers.
- `OUTPUT.selection.tsv` records selected transcript, length basis, protein length and any exclusion reason. A valid GFF3 gene/transcript/CDS hierarchy is required; direct gene-to-CDS or GTF inputs fail explicitly.

Free-text annotation descriptions may contain legacy encoding characters, tabs or semicolons; these are tolerated without using them for feature relationships. Invalid output identifiers are rejected.

## Orthogroups

Run OrthoFinder 2.5.4 with default settings on the exported per-assembly protein FASTAs. Presence across 126 assemblies defines core (101–126), accessory (2–100), and genome-specific (1) orthogroups. Unassigned proteins form singleton groups.


Tests: `python orthogroups/test_extract_proteins.py`.
GFF3 coordinate and phase convention: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
