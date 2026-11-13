# Chromatin accessibility transfer

Training: M. polymorpha control ATAC-seq (GSE290177). Transfer: M. polymorpha heat stress and A. thaliana control/heat-stress DNase-seq (GSE145742).

## Input generation

```bash
python chromatin_accessibility/make_bryocad_atac_lora_input.py \
  --fasta MpTak_v6.fa --narrowpeak Mp_ATAC_Ctrl_peaks.narrowPeak \
  --window_stride 512 --seq_length 512 --write_parquet \
  --outdir results/Mp_ATAC_Ctrl_finetuning_input
```

## Training

Train the ATAC-seq classifier for 4,000 optimization steps:

```bash
python chromatin_accessibility/lora_fine_tune.py train \
  --train_dir=results/Mp_ATAC_Ctrl_finetuning_input/Train.parquet \
  --valid_dir=results/Mp_ATAC_Ctrl_finetuning_input/Validation.parquet \
  --model_name=/path/to/bryocad_small_final \
  --task_type=classification \
  --output_dir=results/atac_fine_tuned_lora_checkpoint_Mp_Ctrl \
  --train_batch_size=8 --eval_batch_size=32 \
  --gradient_accumulation_steps=16 --max_steps=4000 \
  --eval_steps=200 --save_steps=200 --bf16=True
```

The per-device training batch size is 8 with 16 gradient-accumulation steps. BF16 is enabled. Use the final checkpoint at step 4,000 for prediction. An HPC job example is provided in `hpc_reference/finetune_atac.sh`.

## Prediction and evaluation

`hpc_reference/predict_atac.sh` loads checkpoint-4000 and the same local BryoCAD backbone, iterating the validation Parquet files. Match the exact final dataset list before recreating S9. AUROC uses positive-class probabilities; F1 uses argmax without optimized thresholds.

M. polymorpha chr7/chr8 were excluded from training and used for Ctrl validation and HS testing. This is condition transfer on chromosomes absent from training, not evidence that HS loci were absent from validation. Local merged/called peak files and GSM mapping remain source-data inputs to archive.

Program files and adapted job examples are listed with checksums in `SOURCE_FILES.tsv` at the repository root.

## Released adapter

The trained adapter is available at [bryocad/BryoCAD-Small-ATAC-Mpolymorpha-Ctrl-LoRA](https://huggingface.co/bryocad/BryoCAD-Small-ATAC-Mpolymorpha-Ctrl-LoRA). Pair it with BryoCAD-Small. See the model card for checkpoint files and loading instructions.

## HF adapter inference compatibility

Run `bash scripts/run_atac_check.sh` for a pinned-model loading check. The custom Caduceus classifier does not accept every optional keyword that PEFT forwards (for example `attention_mask=None`). The check applies `scripts/peft_forward_compat.py` before constructing the PEFT model, mirroring the forward-wrapper approach in the original fine-tuning pipeline. Only unsupported arguments whose value is `None` are omitted; actual masks or other unsupported values raise an explicit error. These examples use fixed-length unpadded inputs.

A newly initialized `score.weight` message when loading the masked-language-model backbone is expected: the trained score head is restored by the subsequent adapter-loading step. The check compares every saved adapter tensor, including the head, against the loaded state before inference; it fails if any are absent or unequal. The successful output JSON records this verification.

## Readable TSV prediction

`lora_fine_tune.py predict --data_dir input.tsv --output_file predictions.tsv` accepts raw 512-bp sequences, optional 0/1 labels and metadata. It tokenizes in memory, preserves the input columns and adds positive-class probabilities and argmax labels. Existing tokenized Parquet input remains supported. See [TSV examples](../examples/atac/README.md). `bash scripts/run_atac_check.sh /path/to/input.tsv` runs both adapter-restoration and the original prediction CLI. The bundled default TSV has synthetic placeholder labels, not real positive/negative examples.

## Multiple visible GPUs

The ATAC smoke-test runners expose only the first already-visible GPU, preserving scheduler GPU identifiers. This prevents Trainer from automatically enabling DataParallel for a tiny example. The compatibility wrapper now uses the active replica (`self`) rather than capturing the original model's bound forward method. This corrects cross-device dispatch when models are replicated. Replica-dispatch unit tests pass; actual multi-GPU inference has not been validated.
