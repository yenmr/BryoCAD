# BryoCAD

BryoCAD is a DNA language model adapted from PlantCAD2 through continued pretraining on 17 bryophyte genomes spanning hornworts, liverworts and mosses. This repository provides model-use examples and training, evaluation and figure-analysis code.

Repository: https://github.com/yenmr/BryoCAD

## Models

| Model | Hugging Face |
|---|---|
| BryoCAD-Small | [bryocad/BryoCAD-Small](https://huggingface.co/bryocad/BryoCAD-Small) |
| BryoCAD-Large | [bryocad/BryoCAD-Large](https://huggingface.co/bryocad/BryoCAD-Large) |
| GeneCAD head, four bryophytes | [bryocad/GeneCAD-BryoCAD-Small-4bryo](https://huggingface.co/bryocad/GeneCAD-BryoCAD-Small-4bryo) |
| ATAC LoRA, M. polymorpha control | [BryoCAD-Small-ATAC-Mpolymorpha-Ctrl-LoRA](https://huggingface.co/bryocad/BryoCAD-Small-ATAC-Mpolymorpha-Ctrl-LoRA) |

The GeneCAD head is paired with BryoCAD-Small. Use the corresponding model cards for model details. Public model files and their verified revisions are recorded in [model_revisions.json](environments/model_revisions.json).

## Getting started

- [Unified conda installation: BryoCAD + GeneCAD](environments/README.md)
- [Examples and test commands](examples/README.md)
- [Continued pretraining](pretraining/README.md)
- [Protein reconstruction and orthogroups](orthogroups/README.md)
- [Gene-boundary motif recovery](motif_recovery/README.md)
- [Stop-codon prediction analysis](motif_recovery/tts/README.md)
- [Natural and synthetic variants](variant_analysis/README.md)
- [Chromatin accessibility](chromatin_accessibility/README.md)
- [Gene annotation](gene_annotation/README.md)
- [Figure workflows](figures/README.md)

## Evaluation scope

Motif benchmarks include seven model organisms, a 15-species extension (22 total), a 126-assembly bryophyte analysis and a 132-angiosperm TTS extension. Natural SNP analyses use 135 Marchantia polymorpha accessions. Gene annotation Figure 5D evaluates selected whole chromosomes/contigs from 123 bryophyte pangenome assemblies, targeting approximately 50 Mb per assembly. Figures 5A and 5D evaluate CDS features and coding-transcript structure without UTRs; see [evaluation scope](gene_annotation/EVALUATION_SCOPE.md).

Exact recovery concerns gene-boundary sequences, not TF-binding motif discovery. Variant compatibility scores are not experimentally validated deleteriousness labels. Reference annotation quality and representative-isoform selection affect gene-annotation evaluation.

## Reproducibility and availability

Scripts, selected summary tables, four real-sequence motif examples and a synthetic motif input are provided. Additional toy datasets and figure-source inputs are being prepared. GPU training/inference have not been rerun as part of this repository preparation. Some historical workflows require external data and their original software environments; known scorer limitations are documented with the relevant analysis.

Included source files and their checksums are listed in [SOURCE_FILES.tsv](SOURCE_FILES.tsv). HPC job examples contain paths to adapt locally. See [THIRD_PARTY.md](THIRD_PARTY.md) for upstream tools and attribution. BryoCAD code is distributed under [Apache-2.0](LICENSE); bundled upstream components retain their notices. Citation metadata are provided in [CITATION.cff](CITATION.cff).

## Download predicted annotations

[GeneCAD annotations for five bryophyte genomes](https://predict.southerngenomics.org/bryocad/annotation/) are available as GFF files. These are predicted gene annotations, separate from the genome FASTA assemblies and model checkpoints.

## Release and validation

Version **0.1.0**: [release notes](RELEASE_NOTES.md). The unified conda environment and all supplied inference examples passed on RTX PRO 6000 Blackwell, including complete GeneCAD protein refinement. See [recorded outputs](examples/expected/README.md). Other target GPU architectures have not yet been tested.
