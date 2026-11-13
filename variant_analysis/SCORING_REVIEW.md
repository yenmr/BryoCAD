# SNP scoring behavior and regression tests

VCF scoring masks the variant position and uses log(P(ALT)/P(REF)). For a 512-bp context the target index is 255 (zero-based). Output contains scored SNV records; non-SNV records are omitted. Missing flanks are padded with N. The INFO field remains `plantCAD_zero_shot` for compatibility with downstream analysis.

## Input-handling fixes

- Stop indexing the score-record list after it is exhausted, including a trailing non-SNV record.
- Write a header-only VCF for inputs with no SNVs, without loading a model.
- Pad the left and right flanks independently so a short contig still places the variant at the requested target index.
- Declare the score INFO field before creating the VCF writer.

These changes do not alter the model, masking rule, or LLR formula. Regression tests establish identical extracted windows for all positions in a representative long contig and unchanged arithmetic for supplied probabilities. Historical research scores were not recomputed. Full GPU/tokenizer equivalence has not been tested in this release environment.

```bash
python variant_analysis/test_snp_scoring.py
# With PyVCF3 installed:
python variant_analysis/test_snp_vcf_io.py
```

## Input requirements

VCF records must be grouped by chromosome, chromosome names must match FASTA, and REF alleles must agree with the reference genome. REF agreement is not automatically checked. Use valid A/C/G/T SNVs. The released tokenizer must preserve one token per nucleotide and alignment of the masked position. Nonfinite scores can occur if model probabilities underflow and should be inspected before downstream analysis.

Dependencies include PyTorch, transformers, pandas, NumPy, Biopython, tqdm and a package exposing the PyVCF-compatible `vcf` API. The dependency-free regression tests isolate context extraction, score writing and the empty-input path; they do not replace a full model inference test.

VCF serialization and re-reading were tested with PyVCF3 1.0.4.
