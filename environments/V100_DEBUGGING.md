# V100 inference checks

The reported crash occurred after model loading, during the first forward pass:
`LLVM ERROR: Failed to compute parent layout for slice layout.`

This indicates a native compiler/kernel failure, not a failed HF download. The `torch_dtype` deprecation warning is not the fatal error. The previous motif runner selected FP16 automatically on V100; hardware FP16 support does not establish compatibility with every Mamba/Triton kernel. The checks now use FP32 explicitly. Pre-Ampere automatic motif precision also defaults to FP32; explicit `--dtype float16` remains available for controlled comparison.

This is a diagnostic compatibility change, not a confirmed fix on the target V100. The reported stack does not identify the failing kernel or establish whether FP16, Triton version or another kernel configuration is responsible. No environment upgrade/downgrade is performed.

```bash
bash scripts/run_zero_shot_check.sh
bash scripts/run_atac_check.sh
bash scripts/run_genecad_check.sh /path/to/patched-genecad
```

Or run all requested tests with independent failure reporting:

```bash
GENECAD_PYTHON=/path/to/genecad-env/bin/python \
  bash scripts/run_gpu_checks.sh /path/to/patched-genecad
```

The GeneCAD runner creates a temporary copy of the pinned `predict.sh` in its checkout directory, changes the hard-coded BF16 setting to FP32, and removes the temporary file on exit. Original GeneCAD files are unchanged. That checkout directory must be writable. `GENECAD_DTYPE` overrides this test precision; do not select BF16 on V100. Other implementation-specific precision paths are not claimed to have been validated.

All runners enable Python faulthandler and record Triton version, CUDA runtime, GPU compute capability and package versions. A successful import check does not establish kernel compilation compatibility. Downloaded snapshots are reused from the HF cache.

If FP32 still aborts, compare a single donor example in two independent processes:

```bash
bash scripts/run_motif_dtype_probe.sh
```

Return `results/motif_dtype_probe/` and the failing task directory. This records FP32 and FP16 logs and exit codes even when a process aborts. It does not change model weights or training data, clear caches, or silently replace fused kernels. Numerical equivalence between precision modes is not assumed.

Related upstream report of the same error text: https://github.com/state-spaces/mamba/issues/566 (not proof of the cause of this run).
