"""Read human-readable, fixed-length ATAC classification input."""
import csv

def read_atac_tsv(path, sequence_length=512):
    with open(path, newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames or 'sequence' not in reader.fieldnames:
            raise ValueError('TSV requires a sequence column')
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Duplicate TSV column names')
        rows = list(reader)
    if not rows:
        raise ValueError('TSV has no input rows')
    for index, row in enumerate(rows, 2):
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f'Row {index}: inconsistent TSV field count')
        seq = row['sequence'].upper()
        if len(seq) != sequence_length or set(seq) - set('ACGTN'):
            raise ValueError(f'Row {index}: sequence must contain {sequence_length} A/C/G/T/N bases')
        row['sequence'] = seq
        if row.get('label', '') not in ('', '0', '1'):
            raise ValueError(f'Row {index}: optional label must be 0 or 1')
    return rows
