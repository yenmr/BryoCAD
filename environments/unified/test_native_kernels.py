"""Exercise native kernels on the selected GPU before model-level smoke tests."""
import json
import torch
from flash_attn import flash_attn_func
from causal_conv1d import causal_conv1d_fn
assert torch.cuda.is_available(),'Run on an allocated NVIDIA GPU node'
capability=torch.cuda.get_device_capability(0)
assert capability[0]>=8,'This unified profile targets Ampere or newer; not V100'
q=torch.randn(2,128,8,64,device='cuda',dtype=torch.bfloat16)
assert torch.isfinite(flash_attn_func(q,q,q,dropout_p=0.0,causal=False)).all()
x=torch.randn(2,256,128,device='cuda',dtype=torch.bfloat16)
w=torch.randn(256,4,device='cuda',dtype=torch.bfloat16)
assert torch.isfinite(causal_conv1d_fn(x,w)).all()
torch.cuda.synchronize()
print(json.dumps({'gpu':torch.cuda.get_device_name(0),'capability':capability,'flash_attention':'passed','causal_conv1d':'passed','mamba':'validated by subsequent actual BryoCAD model inference'},indent=2))
