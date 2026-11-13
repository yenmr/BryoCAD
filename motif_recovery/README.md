# Gene-boundary motif recovery

Tasks: donor, acceptor, TIS and TTS. Inputs are 8,192-bp oriented windows. The primary metric is exact motif recovery. All target positions are masked simultaneously for BryoCAD/PlantCAD2; A/C/G/T argmax is decoded independently in one pass.

## Generate inputs

Run all commands from the repository root. Create a tab-separated manifest with one row per genome and these three required columns:

- `short_name`: a unique identifier used as the output subdirectory name; do not include `/` or `\`.
- `fa_path`: path to the genome FASTA file, optionally gzip-compressed.
- `gff_path`: path to the corresponding GFF3 annotation, optionally gzip-compressed.

Example (`manifest.tsv`; columns are separated by tabs):

```text
short_name	fa_path	gff_path
M_polymorpha	genomes/MpTak_v6.fa.gz	genomes/MpTak_v6.1r2.gff.gz
A_thaliana	genomes/Ath.fa	genomes/Ath.gff3
```

Relative input paths are resolved against the **manifest directory**, not the shell's working directory. Absolute paths are also accepted. FASTA sequence identifiers must match the GFF3 sequence identifiers.

The parser also accepts `genome_path` as an alias for `fa_path`, and `output` as an alias for `short_name`. If both forms are present, their values must agree. Other columns are optional. The `division` column is required only when using `--angiosperms-only`; accepted values for this filter are `Rosids`, `Asterids`, `Monocots`, `Eudicots`, `Magnoliids`, and `Basal Angio`.

Check the input files before generating sequences:

```bash
python motif_recovery/make_zs_recovery_inputs_from_manifest.py \
  --manifest manifest.tsv \
  --out results/motif_inputs \
  --check-only
```

Generate the four task inputs:

```bash
python motif_recovery/make_zs_recovery_inputs_from_manifest.py \
  --manifest manifest.tsv \
  --out results/motif_inputs
```

For the first row above, outputs are placed in `results/motif_inputs/M_polymorpha/`, including `donor_recovery.tsv`, `acceptor_recovery.tsv`, `tis_recovery.tsv`, and `tts_recovery.tsv`, together with QC summaries. Use `--species M_polymorpha` to process only that manifest identifier.

## Score

```bash
python motif_recovery/score_plantcad2_motif_recovery.py   --input results/motif_inputs/SPECIES/tis_recovery.tsv   --output results/SPECIES_tis_scored.tsv --task tis   --model bryocad/BryoCAD-Small --batch-size 4 --num-workers 0
```

Replace the model with BryoCAD-Large or the relevant PlantCAD2 repository. See `score_evo2_motif_recovery.py --help` and the HPC reference scripts for Evo2 settings. Evo2 forward/reverse-complement comparison must follow the original aggregation rule; do not choose the better direction separately for each site unless that is the intended benchmark.

## Aggregate

```bash
python motif_recovery/process_recovery.py   --results /path/to/scored_results --samples /path/to/samples_list.tsv   --genomes /path/to/genome_summary.tsv   --orthogroups /path/to/Orthogroups_long_classified.tsv   --outdir results/recovery_processed
```

The provided 126-assembly processed tables use **canonical reference filtering**, retain duplicate sites, and report mapping/QC warnings. These are analysis choices already recorded in the source `run_config.json`; they must not be silently generalized to every benchmark. Seven-organism site-generation QC counts are not necessarily final scoring denominators.

The seven-organism merged CSV retains the legacy `motif_accuracy` label; the author states it denotes exact recovery. Preserve the original table and document this interpretation, rather than silently relabeling historical data.

The plotting script also generates exploratory genomic-feature plots. Only panels selected in the figure manifest should enter the manuscript.

## Further analysis of TTS predictions

[Stop-codon prediction analysis](tts/README.md) examines TAA, TAG, and TGA predictions and their relationship to sequence composition. It uses TTS motif-recovery outputs as a downstream analysis.
