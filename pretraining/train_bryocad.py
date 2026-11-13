import os
import time
import subprocess
import argparse
import torch
import random
from datasets import load_from_disk
from transformers import (
    AutoModelForMaskedLM,
    AutoConfig,
    Trainer,
    TrainingArguments,
    TrainerCallback,
    set_seed
)

# ==========================================
# Configure hardware acceleration before training.
# ==========================================
# Enable TF32 for matrix multiplication and convolution.
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True

def get_gpu_utilization():
    """Return GPU compute utilization (%) reported by nvidia-smi."""
    try:
        # Query GPU utilization as numeric values.
        result = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=utilization.gpu", "--format=csv,noheader,nounits"],
            encoding="utf-8"
        )
        # Use the first GPU reported by nvidia-smi.
        utils = [int(x) for x in result.strip().split('\n')]
        return utils[0] if len(utils) > 0 else 0
    except Exception:
        return -1 # Return -1 if utilization is unavailable.

# ==========================================
# Training dashboard: memory and epoch progress.
# ==========================================
class HardwareLoggerCallback(TrainerCallback):
    """Report epoch, step, loss, memory, GPU utilization, and training speed."""
    def __init__(self):
        self.last_time = time.time()
        self.last_step = 0

    def on_log(self, args, state, control, logs=None, **kwargs):
        if logs is not None and "loss" in logs:
            current_time = time.time()
            current_step = state.global_step
            
            # Calculate seconds per step.
            steps_diff = current_step - self.last_step
            time_diff = current_time - self.last_time
            
            if steps_diff > 0 and self.last_step > 0:
                time_per_step = time_diff / steps_diff
                speed_str = f"{time_per_step:.2f} s/it"
            else:
                speed_str = "Calculating..."
                
            self.last_time = current_time
            self.last_step = current_step

            # Read hardware statistics.
            if torch.cuda.is_available():
                max_mem_gb = torch.cuda.max_memory_allocated() / (1024 ** 3)
                gpu_util = get_gpu_utilization()
                util_str = f"{gpu_util}%" if gpu_util >= 0 else "N/A"
                
                current_epoch = state.epoch if state.epoch is not None else 0.0
                
                # Print training progress.
                print(f"[Epoch {current_epoch:.4f} | Step {current_step:5d}] "
                      f"Loss: {logs['loss']:.4f} | "
                      f"VRAM: {max_mem_gb:.2f} GB | "
                      f"GPU Util: {util_str:>3} | "
                      f"Speed: {speed_str}")
                
                # Reset peak memory statistics.
                torch.cuda.reset_peak_memory_stats()

# ==========================================
# Dynamic random cropping and masked-language-model collator.
# ==========================================
class DynamicRandomCropMLMCollator:
    def __init__(self, max_length=8192, mlm_prob=0.15):
        self.max_length = max_length
        self.mlm_prob = mlm_prob
        
        # PlantCAD2 Token Map
        self.PAD_ID = 0
        self.MASK_ID = 1
        self.UNK_ID = 2 # N
        self.VOCAB_DNA = [3, 4, 5, 6] # A, C, G, T

    def __call__(self, features):
        batch_inputs = []
        
        # Randomly crop input sequences.
        for feature in features:
            seq = feature["input_ids"]
            seq_len = len(seq)
            
            if seq_len > self.max_length:
                start_idx = random.randint(0, seq_len - self.max_length)
                crop_seq = seq[start_idx : start_idx + self.max_length]
            else:
                pad_len = self.max_length - seq_len
                crop_seq = seq + [self.PAD_ID] * pad_len
                
            batch_inputs.append(torch.tensor(crop_seq, dtype=torch.long))

        inputs = torch.stack(batch_inputs)
        labels = inputs.clone()

        # Select masked positions with probability 0.15.
        probability_matrix = torch.full(labels.shape, self.mlm_prob)
        
        # Exclude PAD (0) and UNK (2) from masking.
        special_tokens_mask = (inputs == self.PAD_ID) | (inputs == self.UNK_ID)
        probability_matrix.masked_fill_(special_tokens_mask, value=0.0)

        masked_indices = torch.bernoulli(probability_matrix).bool()
        labels[~masked_indices] = -100 # Ignore unselected positions in the loss.

        # Apply the 80/10/10 replacement rule.
        # Replace 80% of selected positions with MASK.
        indices_replaced = torch.bernoulli(torch.full(labels.shape, 0.8)).bool() & masked_indices
        inputs[indices_replaced] = self.MASK_ID

        # Replace 10% with random nucleotides.
        indices_random = torch.bernoulli(torch.full(labels.shape, 0.5)).bool() & masked_indices & ~indices_replaced
        random_dna = torch.randint(3, 7, labels.shape, dtype=torch.long)
        inputs[indices_random] = random_dna[indices_random]

        return {"input_ids": inputs, "labels": labels}

# ==========================================
# Training entry point.
# ==========================================
def parse_args():
    parser = argparse.ArgumentParser(description="BryoCAD Foundation Model Training")
    parser.add_argument("--model_name_or_path", type=str, required=True, help="PlantCAD2 checkpoint path or HF Hub ID")
    parser.add_argument("--data_dir", type=str, required=True, help="Path to local Arrow dataset (e.g., /tmp/...)")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save checkpoints (NAS)")
    parser.add_argument("--per_device_train_batch_size", type=int, default=4)
    parser.add_argument("--gradient_accumulation_steps", type=int, default=8)
    parser.add_argument("--learning_rate", type=float, default=5e-5)
    parser.add_argument("--max_steps", type=int, default=50000)
    parser.add_argument("--save_steps", type=int, default=2000)
    parser.add_argument("--logging_steps", type=int, default=50)
    return parser.parse_args()

def main():
    args = parse_args()
    set_seed(42)

    # Load datasets.
    print(f"\n>>> Loading Arrow Dataset from local NVMe: {args.data_dir}")
    train_dataset = load_from_disk(os.path.join(args.data_dir, "train"))
    val_dataset = load_from_disk(os.path.join(args.data_dir, "val"))
    MAX_EVAL_SAMPLES = 5000
    if len(val_dataset) > MAX_EVAL_SAMPLES:
        # Evaluate at most the first 5,000 validation examples.
        val_dataset = val_dataset.select(range(MAX_EVAL_SAMPLES))
    print(f"Train samples: {len(train_dataset)} | Val samples: {len(val_dataset)}")

    # Load the model with bfloat16 and trust_remote_code.
    print(f"\n>>> Loading Foundation Model: {args.model_name_or_path}")
    config = AutoConfig.from_pretrained(args.model_name_or_path, trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained(
        args.model_name_or_path, 
        config=config, 
        trust_remote_code=True,
        torch_dtype=torch.bfloat16
    )

    # Configure training arguments.
    training_args = TrainingArguments(
        output_dir=args.output_dir,
        overwrite_output_dir=False,
        max_steps=args.max_steps,
        
        # Checkpoint and evaluation schedules.
        eval_strategy="steps",
        eval_steps=args.save_steps,
        save_strategy="steps",
        save_steps=args.save_steps,
        logging_strategy="steps",
        logging_steps=args.logging_steps,
        disable_tqdm=True,

        # Hardware and batch settings.
        per_device_train_batch_size=args.per_device_train_batch_size,
        per_device_eval_batch_size=args.per_device_train_batch_size,
        gradient_accumulation_steps=args.gradient_accumulation_steps,
        bf16=True, 
        dataloader_num_workers=4, # Use CPU workers for sequence cropping.
        
        # Optimizer settings.
        learning_rate=args.learning_rate,
        weight_decay=0.01,
        warmup_steps=1000,
        
        # Retain the best model.
        load_best_model_at_end=True,
        metric_for_best_model="eval_loss",
        save_total_limit=5, 
        
        save_safetensors=False,

        report_to="none"
    )

    # Initialize the Trainer.
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
        data_collator=DynamicRandomCropMLMCollator(max_length=8192, mlm_prob=0.15),
        callbacks=[HardwareLoggerCallback()]
    )

    # Start training.
    print("Starting continued pretraining...")
    try:
        # Resume from the latest checkpoint in the output directory.
        trainer.train(resume_from_checkpoint=True)
    except ValueError:
        print("No checkpoint detected; starting a new training run.")
        trainer.train()

    # Save the final model.
    final_path = os.path.join(args.output_dir, "final_model")
    print(f"\n>>> Training complete. Final model saved to: {final_path}")
    trainer.save_model(final_path)

if __name__ == "__main__":
    main()
