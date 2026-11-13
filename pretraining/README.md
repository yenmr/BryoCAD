# Continued pretraining

Prepare genome FASTA and GFF3 annotations for the [17 selected training genomes](../manifests/training_genomes_17.tsv). Both plain-text and gzip-compressed inputs are supported. GFF3 files must contain `gene` records with unique `ID` attributes and sequence identifiers matching the FASTA headers.

## Build the dataset

Create a tab-separated input manifest with these columns:

```text
species	fasta	gff
M_polymorpha	/path/to/genome.fa.gz	/path/to/annotation.gff3.gz
```

Use one row per selected assembly. Relative input paths are resolved against the manifest directory. The species information table linked above documents genome selection; create this input manifest with the paths to your downloaded files.

From the repository root:

```bash
python pretraining/build_dataset.py \
  --manifest /path/to/pretraining_inputs.tsv \
  --output-dir results/pretraining_dataset
```

The builder retains contigs of at least 20 kb and extracts gene regions with 5 kb flanks on each side. It excludes truncated regions, regions longer than 50 kb, and regions with more than 10% unknown bases. Non-ACGT symbols are treated as unknown; negative-strand regions are reverse-complemented. The output uses PlantCAD2 nucleotide token IDs.

Whole contigs are shuffled globally with seed 42 and assigned to validation until approximately 20% of retained sequence length is reached; remaining contigs form the training partition. Species representation in both partitions is not guaranteed. The output contains `train/`, `val/`, `split_manifest.json`, and `species_map.json`. Use a new output directory for each build.

## Train

```bash
python pretraining/train_bryocad.py \
  --model_name_or_path /path/to/starting_PlantCAD2_checkpoint \
  --data_dir results/pretraining_dataset \
  --output_dir results/pretraining \
  --per_device_train_batch_size 4 \
  --gradient_accumulation_steps 8 \
  --learning_rate 5e-5 \
  --max_steps 50000 \
  --save_steps 2000
```

Choose the Small or Large starting backbone as appropriate. The training script evaluates at most the first 5,000 validation examples. The command specifies a new training run; its step limit is not a statement about the training duration of the released model checkpoints.
