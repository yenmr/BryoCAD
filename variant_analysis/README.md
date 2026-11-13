# Variant analyses

Run all commands from the BryoCAD repository root. Script paths such as `variant_analysis/zero_shot_score.py` are relative to that directory. Replace example input filenames with your data paths; absolute paths are needed only if your data are stored elsewhere.

`zero_shot_score.py` includes input-handling fixes for empty/non-SNV inputs and short contigs. The masked-site scoring formula is unchanged. See [scoring notes](SCORING_REVIEW.md) for input requirements and regression tests.

## Natural SNPs (Figure 4C)

The dataset is M. polymorpha, 135 accessions: https://marchantia.info/variants/data/Marpol_135_accessions/ (PMID 39962240).

`zero_shot_score.py -input-vcf` performs masked-site scoring. `annotate_all_scored_snp_context.py` reads an already-scored VCF and adds genomic context and allele-frequency information; it does **not** calculate model LLRs.

```bash
python variant_analysis/zero_shot_score.py \
  -input-vcf Marpol_135_accessions.vcf.gz -input-fasta MpTak_v6.fa \
  -model bryocad/BryoCAD-Small -output Marpol_135_accessions.BryoCAD_small.scored.vcf \
  -device cuda:0 -batchSize 64 -contextSize 512

python variant_analysis/annotate_all_scored_snp_context.py \
  --vcf Marpol_135_accessions.BryoCAD_small.scored.vcf \
  --fasta MpTak_v6.fa --gff MpTak_v6.1r2.gff.gz \
  --og_table Orthogroups_long_classified.tsv \
  --out Marpol_135_accessions.BryoCAD_small.scored.context.tsv
```

For a model comparison, repeat scoring and context annotation with PlantCAD2, using distinct output filenames. The plot script expects both resulting `.tsv.gz` files in one input directory:

```bash
Rscript variant_analysis/plot_natural_snp_comparison.R /path/to/context_tables results/natural_snp
```

It writes both cumulative top subsets and disjoint rank bands. The visible Figure 4C alone does not establish which output was selected; retain this distinction. Legacy `deleterious_score` fields denote negative LLR, not validated deleterious variants.

## Synthetic SNVs (Figure 4A)

Keep the working directory at the repository root. The commands below assume the input data are in that directory; otherwise provide their paths in the arguments.

```bash
python variant_analysis/generate_core_gene_snvs.py \
  --fasta MpTak_v6.fa --gff MpTak_v6.1r2.gff.gz \
  --og_table Orthogroups_long_classified.tsv --og_species M_polymorpha_v61 \
  --out_prefix core_gene_snvs --max_per_category 10000

python variant_analysis/zero_shot_score.py \
  -input-vcf core_gene_snvs.vcf.gz -input-fasta MpTak_v6.fa \
  -model bryocad/BryoCAD-Small -output core_gene_snvs.scored.vcf \
  -device cuda:0 -batchSize 64 -contextSize 512

Rscript variant_analysis/plot_synthetic_snv_scores.R core_gene_snvs.annotation.tsv core_gene_snvs.scored.vcf
```

LLR is log P(alt) - log P(ref), at the masked position. The supplied generation command's local backbone was `BryoCAD_small10`; verify identity with the HF release before claiming identical scores.

## Synthetic deletions (Figure 4B; MLM only)

```bash
python variant_analysis/make_synthetic_deletion_input.py \
  --manifest genome_gff_manifest.tsv --out synthetic_deletion_input.tsv \
  --seq_length 8192 --tis_upstream_del_len 100 --cds_inframe_del_len 9 \
  --cds_frameshift_del_len 10 --cds_flank 10 --max_per_type_per_genome 1000 \
  --max_n_ratio 0.05 --seed 42

python variant_analysis/score_synthetic_deletion_sv_mlm.py \
  --input synthetic_deletion_input_core.tsv --output synthetic_deletion_core_scored_mlm.tsv \
  --model bryocad/BryoCAD-Small --device cuda:0 --batch_size 8 \
  --seq_length 8192 --flank_size 10
```

Use the actual core-gene subset from the study. The 12,000 events/7,881 loci describe the generated dataset, not Supplementary Table S1. MLM deletion scores are based on retained flanking positions and are not the single-position SNP formula. Unmasked deletion scoring is excluded.
