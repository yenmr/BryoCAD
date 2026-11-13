# SNP scoring

## Inputs

[input.vcf](input.vcf) contains 20 biallelic SNPs with `FILTER=PASS`. The [reference FASTA](../gene_annotation/input.fa.gz) covers *Marchantia polymorpha* chr2:1–200000. FASTA and VCF sequence identifiers match, and coordinates require no offset.

Each SNP has at least 4,096 bp available on both sides within the reference fragment. The 8,192-bp scoring window contains at most 1% N bases. REF alleles match the reference sequence.

## Run

Activate the [conda environment](../../environments/README.md), then run from the repository root:

```bash
bash scripts/run_zero_shot_check.sh
```

The script downloads BryoCAD-Small and runs both motif recovery and SNP scoring. SNP scoring uses an 8,192-bp window, with the masked variant at zero-based index 4,095.

## Output

```text
results/zero_shot_check/snps.vcf
```

The output preserves the 20 input variants and adds `plantCAD_zero_shot` to the INFO field:

```text
LLR = log P(ALT | masked context) - log P(REF | masked context)
```

A positive score indicates that the model assigns higher probability to ALT; a negative score favors REF. These scores alone do not establish a variant's biological effect.

## Use your own variants

Provide a VCF with the eight standard columns and a matching reference FASTA. For SNP examples, use single-base A/C/G/T REF and ALT alleles. VCF positions are 1-based; if you extract a reference fragment starting after chromosome position 1, convert variant positions to fragment-relative coordinates.

Accession genotype columns are unnecessary for LLR scoring. Retain them separately if you also intend to calculate allele frequencies.
