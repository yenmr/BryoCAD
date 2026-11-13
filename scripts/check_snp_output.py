import argparse, math
import vcf
p=argparse.ArgumentParser(); p.add_argument('vcf'); p.add_argument('--input', required=True); a=p.parse_args()
with open(a.input) as f: expected=list(vcf.Reader(f))
with open(a.vcf) as f: rows=list(vcf.Reader(f))
key=lambda r:(r.CHROM,r.POS,r.REF,tuple(str(x) for x in r.ALT))
assert rows and list(map(key,rows))==list(map(key,expected)), 'Output variants differ from input'
for r in rows:
 values=r.INFO['plantCAD_zero_shot']
 if not isinstance(values,list): values=[values]
 assert len(values)==len(r.ALT) and all(math.isfinite(float(x)) for x in values)
print(f'{len(rows)} finite SNP scores validated')
