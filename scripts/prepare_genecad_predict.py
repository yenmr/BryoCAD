"""Set explicit precision in a temporary GeneCAD wrapper without editing upstream."""
import argparse
import re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('output',type=Path)
p.add_argument('--dtype',choices=['float32','float16','bfloat16'],default='float32');a=p.parse_args()
s=a.source.read_text();old='DTYPE="bfloat16"'
if s.count(old)!=1:raise SystemExit('Expected one DTYPE="bfloat16" assignment in pinned predict.sh; inspect this upstream version before adapting.')
s=s.replace(old,'DTYPE="'+a.dtype+'"')
s,count=re.subn(r'(local\s+max_attempts=)\d+', r'\g<1>1',s)
if count!=1:raise SystemExit('Expected one max_attempts setting in pinned predict.sh; inspect this upstream version.')
a.output.write_text(s)
print('GeneCAD wrapper dtype:',a.dtype)
