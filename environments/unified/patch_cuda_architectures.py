"""Extend the user-validated Mamba/causal-conv1d setup patch to target GPUs."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('files',nargs='+',type=Path);p.add_argument('--report',type=Path,required=True);a=p.parse_args()
marker='    # HACK: The compiler flag -D_GLIBCXX_USE_CXX11_ABI'
records=[]
for path in a.files:
    original=path.read_text();before=hashlib.sha256(original.encode()).hexdigest()
    additions=[arch for arch in ['80','89','90','120'] if f'arch=compute_{arch},code=sm_{arch}' not in original]
    if additions:
        if original.count(marker)!=1:raise SystemExit(f'Unexpected upstream setup.py; no unique insertion point: {path}')
        snippet='\n    # BryoCAD CUDA 12.8 target architectures; preserve upstream flags.\n    if bare_metal_version >= Version("12.8"):\n'
        for arch in additions:snippet+=f'        cc_flag.extend(["-gencode", "arch=compute_{arch},code=sm_{arch}"])\n'
        path.write_text(original.replace(marker,snippet+'\n'+marker,1))
    records.append({'file':str(path),'before_sha256':before,'after_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'added_sm':additions})
a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
