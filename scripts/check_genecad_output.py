import argparse,gzip,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--fasta',type=Path,required=True);a=p.parse_args()
lengths={}
with (gzip.open(a.fasta,'rt') if a.fasta.suffix=='.gz' else a.fasta.open()) as f:
 for line in f:
  if line.startswith('>'): name=line[1:].split()[0];lengths[name]=0
  else: lengths[name]+=len(line.strip())
files=[x for x in a.directory.rglob('*') if x.is_file() and x.name.endswith(('.gff','.gff3','.gff.gz','.gff3.gz'))]
assert files,'No GFF output was produced'
summary={}
for path in files:
 count=0
 with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:
  for line in f:
   if line.startswith('#') or not line.strip():continue
   fields=line.rstrip().split('\t');assert len(fields)==9,path
   assert fields[0] in lengths and 1<=int(fields[3])<=int(fields[4])<=lengths[fields[0]],(path,fields[:5])
   count+=1
 summary[str(path.relative_to(a.directory))]=count
(a.directory.parent/'gff_check.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
final=[n for n in summary if n.endswith('_GeneCAD_final.gff')]
assert final and all(summary[n]>0 for n in final), 'Expected a nonempty final GFF for the complete-gene example'
print('Nonempty final GFF and coordinate checks passed; this is not an accuracy benchmark.')
