# GeneCAD version and supplied modifications

Upstream: https://github.com/plantcad/genecad
Release: v0.4.0; locally verified commit `1a9baf84d61152c20a0cd910694316cd88d1d374`.

The author supplied the actual `train.sh` and `src/modeling.py`. Their complete differences from the pinned upstream files are now preserved in `upstream_modifications/genecad-v0.4.0-bryocad.patch`. The original files and upstream Apache-2.0 license are also included. No further full-checkout request is needed for these known modifications.

```bash
git clone https://github.com/plantcad/genecad.git
cd genecad
git checkout 1a9baf84d61152c20a0cd910694316cd88d1d374
git apply /path/to/BryoCAD/gene_annotation/upstream_modifications/genecad-v0.4.0-bryocad.patch
```

The patch was checked, applied to clean copies of the pinned files, and verified byte-for-byte against both supplied files. This checks patch reconstruction, not a new GPU training run.

## Actual changes and effective wrapper settings

| Setting | Supplied internal code | Effective 4bryo wrapper |
|---|---|---|
| Batch size / effective batch | Defaults 8 / 40 | 8 / 40 |
| Learning rate | Default 1e-4 | 1e-4 |
| Token embedding | 128 for plant | 128 |
| Base frozen | yes for plant | explicitly yes |
| Automatic class weights | default no | explicitly **yes**, overriding default |
| Base model | default PlantCAD2 legacy Small | explicitly local BryoCAD-Small |
| Epochs | plant default 1 | explicitly 5 |
| Scheduler | none (constant) | same |
| Weight decay | 1e-5 | same |
| Validation interval | 1.0 | same |
| Train-eval frequency | 100 | same |

The run-name default also changes to a legacy-head label; the 4bryo wrapper retains an older 2bryo run-name label, although it selects four species. Preserve the actual species arguments as evidence.

Published head: https://huggingface.co/bryocad/GeneCAD-BryoCAD-Small-4bryo . Immutable revisions/local checkpoint mapping can be recorded at final release. `capture_upstream_state.sh` is optional archival tooling, not a current missing-file requirement.
