import argparse,csv,math
p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output');a=p.parse_args()
def read(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
source,output=read(a.input),read(a.output)
assert source and len(source)==len(output)
for x,y in zip(source,output):
    for key in x:assert y[key]==(x[key].upper() if key=='sequence' else x[key]),key
    probability=float(y['probability_positive'])
    assert math.isfinite(probability) and 0<=probability<=1
    assert y['predicted_label'] in ('0','1')
print(f'PASS: {len(output)} TSV predictions with input metadata preserved')
