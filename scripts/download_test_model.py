"""Download one pinned HF snapshot; print its local path for shell runners."""
import argparse,json
from pathlib import Path
from huggingface_hub import snapshot_download
p=argparse.ArgumentParser();p.add_argument('name');a=p.parse_args()
rows=json.loads((Path(__file__).resolve().parents[1]/'environments/model_revisions.json').read_text())
row=next(x for x in rows if x['repository']=='bryocad/'+a.name)
print(snapshot_download(repo_id=row['repository'],revision=row['revision']))
