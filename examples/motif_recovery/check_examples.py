"""Check toy inputs and, optionally, output labels without model dependencies."""
import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--results', type=Path)
args = parser.parse_args()
base = Path(__file__).resolve().parent
for task in ['donor', 'acceptor', 'tis', 'tts']:
    with (base / (task + '.tsv')).open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert len(rows) == 4, task
    assert {r['direction'] for r in rows} == {'+', '-'}, task
    for row in rows:
        sequence = row['sequence']; start = int(row['motif_start']); end = int(row['motif_end'])
        assert len(sequence) == 8192
        assert sequence[start:end] == row['motif']
        assert end - start == (2 if task in ['donor', 'acceptor'] else 3)
    if args.results:
        with (args.results / (task + '.tsv')).open() as f:
            predicted = list(csv.DictReader(f, delimiter='\t'))
        assert len(predicted) == len(rows)
        for source, result in zip(rows, predicted):
            assert source['geneid'] == result['geneid']
            assert source['motif'] == result['motif']
            assert result['exact_recovery'].lower() in ['true','false']
            assert (result['exact_recovery'].lower() == 'true') == (result['predicted_motif'] == source['motif'])
    print(task + ': 4 records validated')
