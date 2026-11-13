"""Load the published adapter and trained head and run two synthetic inputs."""
import argparse,json
from pathlib import Path
import torch
from transformers import AutoTokenizer,AutoModelForSequenceClassification
from peft import PeftModel, get_peft_model_state_dict
from peft_forward_compat import adapt_forward_for_peft
from safetensors import safe_open
p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--adapter',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
with safe_open(str(Path(a.adapter)/'adapter_model.safetensors'),framework='pt',device='cpu') as f:
 saved_weights={k:f.get_tensor(k) for k in f.keys()}
 head_keys=[k for k in saved_weights if '.score.' in k or '.classifier.' in k]
assert head_keys, 'No trained classification-head weights found in adapter'
tok=AutoTokenizer.from_pretrained(a.base,trust_remote_code=True)
base=AutoModelForSequenceClassification.from_pretrained(a.base,trust_remote_code=True,num_labels=2,torch_dtype=torch.float32)
# The original LoRA pipeline also wraps the custom forward before PEFT loading.
adapt_forward_for_peft(base)
model=PeftModel.from_pretrained(base,a.adapter)
loaded_weights=get_peft_model_state_dict(model)
for key, expected in saved_weights.items():
    if key not in loaded_weights:
        raise RuntimeError('Saved adapter tensor was not restored: '+key)
    actual=loaded_weights[key].detach().cpu().float()
    if not torch.equal(actual, expected.float()):
        raise RuntimeError('Restored adapter tensor differs from saved tensor: '+key)
print(f'Verified {len(saved_weights)} restored adapter tensors, including task head: {head_keys}', flush=True)
model=model.to('cuda').eval()
sequences=['ACGT'*128,'ATGC'*128]
inputs=tok(sequences,return_tensors='pt',add_special_tokens=False,return_attention_mask=False,return_token_type_ids=False)['input_ids'].to('cuda')
assert tuple(inputs.shape)==(2,512)
with torch.inference_mode():logits=model(input_ids=inputs).logits
assert tuple(logits.shape)==(2,2) and torch.isfinite(logits).all()
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps({'head_keys':head_keys,'verified_adapter_tensors':len(saved_weights),'saved_weights_restored_exactly':True,'logits':logits.float().cpu().tolist(),'positive_probabilities':logits.float().softmax(-1)[:,1].cpu().tolist(),'purpose':'synthetic loading/inference check, not biological evaluation'},indent=2))
print('ATAC adapter and task head loaded; two finite predictions written')
