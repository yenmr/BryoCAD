# Fixed-policy annotation location evaluator

All genes, transcripts and predictions are retained. No consistency-based gene
filtering, ignored loci, prediction removal, sequence motifs or BUSCO. Python
standard library only; FASTA is not required. Input GFFs are never modified.

```bash
python ../prediction_pangenome/evaluate_compatible.py \
  --ref Araport11_genecad.gff3.gz \
  --pred bryocad_small_d768_4bryo_3species/Ath/Ath_GeneCAD_final.gff \
  --structure-policy utr \
  --metrics-tsv Ath.utr.metrics.tsv \
  --output Ath.utr.evaluation.txt \
  --audit-json Ath.utr.audit.json
```

## Three fixed geometry rules

`--structure-policy` accepts `utr` (default), `exon`, or `intersection`. The same
rule is applied to reference and prediction, to all organisms and all isoforms.
A rule changes effective feature geometry, never which genes are evaluated.
CDS coordinates are always retained and CDS base/segment/chain metrics are
independent of this option.

When a coding transcript has BOTH explicit exon and UTR records:

| Option | Effective UTRs | Effective exons |
|---|---|---|
| utr | Explicit 5-prime/3-prime UTRs | Union of original CDS plus effective UTRs |
| exon | Terminal exonic regions outside the outer CDS boundaries | Explicit exons, preserving any CDS bases outside exons |
| intersection | Explicit UTR intersected with independently exon+CDS-inferred UTR, separately by UTR type | Union of original CDS plus intersected UTRs |

`intersection` is a base-interval intersection, not intersection of gene IDs or
whole exons. It never clips CDS or removes genes. Negative-strand UTR directions
are reversed appropriately. In all policies, overlapping CDS bases take priority
and are removed from effective UTR intervals. Warnings/raw feature statistics
remain available in audit.

For a positive-strand example with CDS 11–90, exon 1–100, explicit 5-prime UTR
1–10 and explicit 3-prime UTR 96–110:

- utr: effective 3-prime UTR 96–110;
- exon: effective 3-prime UTR 91–100;
- intersection: effective 3-prime UTR 96–100.

The three definitions are not claims about which annotation is biologically
correct. A policy should be chosen once and used for every compared model.

If there is no pair of independent exon/UTR representations to reconcile:

- Exon+CDS without explicit UTR: derive terminal UTRs from exon+CDS for all options.
- CDS+UTR without explicit exon: reconstruct exons from CDS+UTR for all options;
  this does not independently validate UTRs.
- Reference with CDS only and no exon: coding metrics are available; full-exon/
  UTR comparison metrics are NA. No UTR is invented from mRNA bounds.
- Predicted CDS-only models are retained as emitted, with no predicted UTR.
- Transcripts without CDS retain available exon/explicit UTR information;
  coding-based UTR inference is not applicable.

## Retaining prediction information when reference lacks UTR

The TSV has `section metric status reason TP FP FN precision recall f1 ref_count
pred_count`. `pred_count` is retained even for unavailable full-structure/UTR
comparisons; for base_level metrics its unit is bases. If reference has only
CDS and prediction has 500 UTR bases, the UTR row contains pred_count=500,
ref_count=NA, status=NA and TP/FP/FN/P/R/F1=NA. Those bases are not incorrectly
called FP and the prediction is not deleted. Raw feature record counts, raw
feature union bp and full transcript count remain in audit JSON. Existing
matrix columns are retained; ref_count/pred_count are appended.

## Stable matrix groups

- `CDS/cds_ex`: exact CDS segments; `CDS/tx`: entire CDS chains;
  `CDS/locus`: gene-level any-isoform CDS-chain matching.
- `coding_structure/intron`, `ichain`: coding-segment gaps and complete gap chains.
- `full_structure/exon`, `intron`, `ichain`, `locus_ic`, `full_tx`: resolved full
  exon structure metrics. full_tx requires all exon boundaries including ends.
- `base_level/exonic`, `cds`, `five_prime_utr`, `three_prime_utr`: base overlap.
  There is no duplicate full_structure/base or coding_structure/base row.
- `sites/CDS.5prime`, `CDS.3prime`, `Coding.Junc.Donor`, `Coding.Junc.Acceptor`,
  `Full.Junc.Donor`, `Full.Junc.Acceptor`: coordinate-only site matching.
  Duplicate TIS/TTS placeholder rows are removed. CDS.5prime/CDS.3prime are
  the coordinate-only equivalents of the original evaluator TIS/TTS metrics;
  no motif verification occurs.

Base intervals are unioned by (sequence ID, strand), counting overlapping
isoforms once within a region. Reference information availability and undefined
denominators are represented with NA. Exact gene/transcript matching is
one-to-one. CDS phase and stop-codon inclusion conventions are not adjusted.
All-zero reference/prediction region sets have undefined P/R/F1, not perfect F1.

## Audit and compatibility

Audit contains `structure_policy`, resolution-source counts, raw record/bp counts,
input warnings and `gene_selection='none; all genes and predictions retained'`.
No exclusion file is generated. Old excluded-gene lists and filtered matrices are
historical outputs and must not be applied to this version.

Reference-only inspection: `--ref reference.gff3.gz --audit-only --audit-json audit.json`.
`--seqid` deliberately restricts both inputs; otherwise differing sequence-ID sets
fail instead of silently dropping chromosomes. GFF3 gzip, geneID hierarchy links,
shared Parent records and out-of-order features are supported. Invalid/missing
hierarchy, inconsistent chromosome/strand links or invalid coordinates can stop
parsing, but no invalid model is silently selected out. GTF is not supported.
`--fasta` is optional coordinate-bounds checking only. `--skip-busco`/`--cpu` are
accepted for old commands; `--reference-structure` is deprecated and ignored.

## Validation

21 regression tests pass. Tests cover all three policies, negative strands,
unchanged one-source geometry, preservation of all genes, retained prediction
UTR counts with CDS-only reference, interval unions, NA handling and exact matching.
The actual Ath data ran under all three options. Each retained 27,561 reference
transcripts and 22,940 predicted transcripts; CDS-chain confusion counts were
identical across options. Ath matrices are saved under Ath.fixed_utr,
Ath.fixed_exon and Ath.fixed_intersection names.
