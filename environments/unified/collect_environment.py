"""Record packages, conda builds, CUDA compiler/driver and editable-source metadata."""
import argparse,importlib.metadata,json,os,platform,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
report={'python':sys.version,'executable':sys.executable,'prefix':sys.prefix,'platform':platform.platform(),'CUDA_HOME':os.environ.get('CUDA_HOME'),'TORCH_CUDA_ARCH_LIST':os.environ.get('TORCH_CUDA_ARCH_LIST'),'packages':sorted([{'name':d.metadata.get('Name','unknown'),'version':d.version} for d in importlib.metadata.distributions()],key=lambda x:x['name'].lower())}
try:
    import torch
    report['torch_cuda']=torch.version.cuda;report['compiled_architectures']=torch.cuda.get_arch_list()
    if torch.cuda.is_available():report['gpu']=torch.cuda.get_device_name(0);report['capability']=torch.cuda.get_device_capability(0)
except Exception as e:report['torch_error']=str(e)
(a.output/'environment.json').write_text(json.dumps(report,indent=2))
commands={'nvcc.txt':['nvcc','--version'],'nvidia-smi.txt':['nvidia-smi'],'conda-explicit.txt':['conda','list','--explicit','-p',sys.prefix],'conda-list.json':['conda','list','--json','-p',sys.prefix],'pip-check.txt':[sys.executable,'-m','pip','check']}
statuses={}
for name,command in commands.items():
    if not shutil.which(command[0]):
        statuses[name]={'available':False};continue
    result=subprocess.run(command,capture_output=True,text=True)
    (a.output/name).write_text(result.stdout)
    if result.stderr:(a.output/(name+'.stderr.txt')).write_text(result.stderr)
    statuses[name]={'available':True,'exit_code':result.returncode}
(a.output/'command-status.json').write_text(json.dumps(statuses,indent=2))
print('Saved environment inventory to',a.output)
