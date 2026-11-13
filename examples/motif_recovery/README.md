# Gene-boundary motif recovery examples

Four real-sequence records per task (donor, acceptor, TIS and TTS) were selected from the author's 250-site-per-task speed-test samples. Each record contains an 8,192-bp sequence and zero-based, half-open motif coordinates. Both source strands are represented. Sequences are already oriented as expected by the scorer; do not reverse-complement them again.

Run from the repository root in a model-compatible environment:

```bash
bash examples/motif_recovery/run.sh
```

The default model is `bryocad/BryoCAD-Small`. Set `MODEL=/path/to/model` to use a local checkpoint. The script runs all four tasks and writes `results/toy_motif_recovery/{task}.tsv`. Each output should have four rows and include `predicted_motif` and `exact_recovery`. Exact recovery is true only if every predicted base equals the reference motif. Small samples are for demonstrating execution, not estimating performance.

Input validation requires only Python's standard library:

```bash
python examples/motif_recovery/check_examples.py
```

After model inference, validate the output structure and exact-recovery labels:

```bash
python examples/motif_recovery/check_examples.py --results results/toy_motif_recovery
```

`expected_summary.tsv` records expected row counts, sequence lengths and motif widths. `selection_manifest.json` records the source-row indices (one-based data rows, excluding header) and hashes; `source_manifest.json` records the original sampling provenance. Original row metadata are retained. Sequence data retain their source terms.
