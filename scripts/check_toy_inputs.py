"""Validate bundled biological example inputs without external dependencies."""
import csv
import gzip
from pathlib import Path
root = Path(__file__).resolve().parents[1]
with gzip.open(root/'examples/gene_annotation/input.fa.gz', 'rt') as handle:
    lines = handle.read().splitlines()
assert lines[0].startswith('>')
chrom = lines[0][1:].split()[0]
sequence = ''.join(lines[1:]).upper()
assert len(sequence) == 200000
records = [line.split('\t') for line in (root/'examples/snp/input.vcf').read_text().splitlines() if not line.startswith('#')]
assert len(records) == 20
for record in records:
    pos = int(record[1]) - 1
    assert record[0] == chrom and sequence[pos] == record[3]
    assert len(record[3]) == len(record[4]) == 1
    assert record[4] in 'ACGT' and record[4] != record[3]
    assert min(pos, len(sequence)-pos-1) >= 4096
    window = sequence[pos-4095:pos-4095+8192]
    assert len(window) == 8192 and window.count('N') / len(window) <= .01
with gzip.open(root/'examples/gene_annotation/reference.gff3.gz', 'rt') as handle:
    for line in handle:
        if line.startswith('#') or not line.strip():
            continue
        fields = line.rstrip().split('\t')
        assert len(fields) == 9 and fields[0] == chrom
        assert 1 <= int(fields[3]) <= int(fields[4]) <= len(sequence)
with (root/'examples/atac/input_template.tsv').open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
assert len(rows) == 10
for row in rows:
    assert len(row['sequence']) == 512 and set(row['sequence'].upper()) <= set('ACGTN')
    assert row['label'] in ('0', '1')
print('GeneCAD FASTA/GFF, SNP REF/context, and ATAC TSV checks passed')
