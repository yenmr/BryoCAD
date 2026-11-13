# Gene-boundary motif recovery

Tasks: donor, acceptor, TIS and TTS. Inputs are 8,192-bp oriented windows. The primary metric is exact motif recovery. All target positions are masked simultaneously for BryoCAD/PlantCAD2; A/C/G/T argmax is decoded independently in one pass.

## Generate inputs

Resolve `genome_path` and `gff_path` in a copy of the source manifest to your local files. Then:

```bash
python motif_recovery/make_zs_recovery_inputs_from_manifest.py   --manifest /path/to/resolved_manifest.tsv --out results/motif_inputs
```

Use `--check-only` to check input paths first. Original manifests retain historical paths as provenance, not working download URLs.

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
