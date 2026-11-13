"""Make transparent macro summaries from the bundled seven-species table."""
from pathlib import Path
import csv
from collections import defaultdict
ROOT = Path(__file__).resolve().parents[1]
source = ROOT / 'figures/source_data/model_organisms/by_task_and_species_merged.csv'
models = {'Evo2-1B':'Evo2-1B','Evo2-7B':'Evo2-7B','Evo2-20B':'Evo2-20B','PlantCAD2-L':'PlantCAD2-Large','BryoCAD-L':'BryoCAD-Large','plantCAD2_Small':'PlantCAD2-Small','bryoCAD_Small':'BryoCAD-Small'}
with source.open() as f: rows=list(csv.DictReader(f))
out=ROOT/'results/model_organisms'; out.mkdir(parents=True,exist_ok=True)
long=[]; groups=defaultdict(list)
for r in rows:
    for col,name in models.items():
        value=float(r[col])
        if not 0 <= value <= 1: raise ValueError(f'Invalid recovery: {value}')
        long.append(dict(species=r['Species'],task=r['Task'],model=name,exact_recovery=value,source_sites=r['Sites']))
        groups[name,r['Task']].append(value)
        groups[name,'ALL_TASKS'].append(value)
for filename,data in [('exact_recovery_long.tsv',long),('macro_summary.tsv',[dict(model=m,task=t,n_cells=len(v),mean_exact_recovery=sum(v)/len(v)) for (m,t),v in groups.items()])]:
    with (out/filename).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
print(f'Wrote {len(long)} species-task-model values to {out}')
