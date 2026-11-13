"""Record the active environment and check CUDA inference prerequisites."""
import argparse, importlib, importlib.metadata, json, platform, subprocess, sys, os
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--profile',choices=['zero-shot','atac','all'],default='all');p.add_argument('--require-cuda',action='store_true');p.add_argument('--output',type=Path);a=p.parse_args()
report={'platform':platform.platform(),'python':sys.version,'executable':sys.executable,'CUDA_VISIBLE_DEVICES':os.environ.get('CUDA_VISIBLE_DEVICES'),'packages':{},'errors':[]}
for name in ['torch','triton','transformers','mamba-ssm','causal-conv1d','peft','huggingface-hub','PyVCF3','numpy','pandas']:
 try:report['packages'][name]=importlib.metadata.version(name)
 except importlib.metadata.PackageNotFoundError:report['packages'][name]=None
if a.require_cuda:
 modules=['torch','transformers','mamba_ssm','causal_conv1d']
 if a.profile in ['zero-shot','all']:modules.append('vcf')
 if a.profile in ['atac','all']:modules.append('peft')
 for name in modules:
  try:importlib.import_module(name)
  except Exception as e:report['errors'].append(name+': '+str(e))
 try:
  import torch
  from cuda_preflight import check_cuda_runtime
  cuda_report, error = check_cuda_runtime(torch)
  report.update(cuda_report)
  if error:report['errors'].append(error)

 except ImportError:pass
if a.output:
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2))
 # uv environments need not contain pip. Record installed distributions directly.
 packages = sorted(
     [{'name': d.metadata.get('Name', 'unknown'), 'version': d.version}
      for d in importlib.metadata.distributions()], key=lambda d:d['name'].lower())
 a.output.with_suffix('.packages.json').write_text(json.dumps(packages, indent=2))
print(json.dumps(report,indent=2))
if report['errors']:raise SystemExit(1)
