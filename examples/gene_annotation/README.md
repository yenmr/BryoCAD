# Gene annotation

## Inputs

- [input.fa.gz](input.fa.gz): a 200-kb *Marchantia polymorpha* fragment covering chr2:1–200000, in genomic orientation.
- [reference.gff3.gz](reference.gff3.gz): annotations for 13 complete genes and their descendants; genes crossing fragment boundaries are omitted.
- [source.json](source.json): assembly information and coordinate conventions.

Only the FASTA is used for prediction. The reference GFF3 is provided for inspecting the resulting annotations.

## Run

Activate the [conda environment](../../environments/README.md), then run from the repository root:

```bash
export GENECAD_PYTHON="$CONDA_PREFIX/bin/python"
bash scripts/run_genecad_check.sh "$HOME/software/genecad-bryocad"
```

Replace the GeneCAD source path if you installed it elsewhere. The script downloads the BryoCAD-Small backbone and GeneCAD head, reads the compressed FASTA directly, and runs annotation followed by protein refinement. Additional protein models are downloaded as needed.

## Output

The final annotation is written to:

```text
results/genecad_check/predictions/toy_mpolymorpha_GeneCAD_final.gff
```

The same directory contains intermediate raw and ORF-adjusted GFF files. Output coordinates are 1-based and relative to the input FASTA. For this fragment, they also match the original chr2 coordinates because the fragment starts at position 1.

The runner checks that the final GFF is nonempty and its coordinates lie within the input sequence. It archives earlier prediction directories when rerun.
