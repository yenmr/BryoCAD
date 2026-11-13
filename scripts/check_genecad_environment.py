"""Fail early on an incompatible GeneCAD interpreter or missing dependencies."""
import argparse,importlib,importlib.metadata,json,os,platform,re,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--checkout',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
report={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'checkout':str(a.checkout),'packages':{},'errors':[]}
project=a.checkout/'pyproject.toml'
if not project.exists():report['errors'].append('GeneCAD pyproject.toml is missing')
else:
    match=re.search(r'requires-python\s*=\s*["\']([^"\']+)',project.read_text())
    if match:
        report['requires_python']=match.group(1)
        try:
            from packaging.specifiers import SpecifierSet
            if platform.python_version() not in SpecifierSet(match.group(1)):
                report['errors'].append('Selected Python does not satisfy GeneCAD requires-python '+match.group(1))
        except ImportError:report['errors'].append('Missing packaging: unable to check Python requirements')
modules={'torch':'torch','transformers':'transformers','lightning':'lightning','xarray':'xarray','zarr':'zarr','numpy':'numpy','pandas':'pandas','Bio':'biopython','numba':'numba','omegaconf':'omegaconf','sklearn':'scikit-learn','scipy':'scipy','huggingface_hub':'huggingface-hub','mamba_ssm':'mamba-ssm','causal_conv1d':'causal-conv1d'}
for module,distribution in modules.items():
    try:
        importlib.import_module(module)
        report['packages'][distribution]=importlib.metadata.version(distribution)
    except Exception as e:report['errors'].append(module+': '+str(e))
try:
    from zarr.errors import GroupNotFoundError
except Exception as e:report['errors'].append('GeneCAD Zarr API check: '+str(e))
if 'torch' in sys.modules:
    torch=sys.modules['torch']
    from cuda_preflight import check_cuda_runtime
    cuda_report, error = check_cuda_runtime(torch)
    report.update(cuda_report)
    if error:report['errors'].append(error)

a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
if report['errors']:
    print('\nSTOP: GeneCAD environment check failed before downloads/pipeline. Inspect the errors above: select the GeneCAD environment for missing dependencies; for CUDA failure use a GPU/build-compatible environment. Do not change the working BryoCAD environment to repair GeneCAD.',file=sys.stderr)
    raise SystemExit(1)
