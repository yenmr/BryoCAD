"""CPU checks for FASTA/GFF3 pretraining input generation."""
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('builder', Path(__file__).resolve().parents[1] / 'pretraining/build_dataset.py')
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


class InputTests(unittest.TestCase):
    def test_coordinates_orientation_and_filters(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            with gzip.open(p / 'genome.fa.gz', 'wt') as out:
                out.write('>chr1\nAACGTNACGTAA\n>chr2\nAAAAAAAAAAAA\n')
            (p / 'genes.gff3').write_text('chr1\tt\tgene\t3\t4\t.\t+\t.\tID=plus\nchr1\tt\tgene\t3\t4\t.\t-\t.\tID=minus\nchr1\tt\tgene\t1\t2\t.\t+\t.\tID=edge\nchr1\tt\tgene\t6\t6\t.\t+\t.\tID=unknown\n')
            (p / 'inputs.tsv').write_text('species\tfasta\tgff\nspecies1\tgenome.fa.gz\tgenes.gff3\n')
            rows = builder.read_manifest(p / 'inputs.tsv')
            split = builder.plan_split(rows, 1, .2, 42)
            self.assertEqual(split, builder.plan_split(rows, 1, .2, 42))
            self.assertEqual(len(split['train']), 1)
            self.assertEqual(len(split['val']), 1)
            records = list(builder.gene_generator(rows, [{'species': 'species1', 'chrom': 'chr1'}], {'species1': 0}, flank=1, max_n_ratio=0))
            self.assertEqual([r['gene_uid'] for r in records], ['plus', 'minus'])
            self.assertEqual(records[0]['input_ids'], [3, 4, 5, 6])
            self.assertEqual(records[1]['input_ids'], [3, 4, 5, 6])
            self.assertEqual(records[0]['roi_len'], 4)
            self.assertEqual(list(builder.gene_generator(rows, [{'species': 'species1', 'chrom': 'chr1'}], {'species1': 0}, flank=1, max_roi_len=2)), [])

    def test_reverse_complement(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / 'f.fa').write_text('>c\nAACGTA\n')
            (p / 'g.gff').write_text('c\tt\tgene\t2\t4\t.\t-\t.\tID=g\n')
            rows = [{'species': 's', 'fasta': str(p/'f.fa'), 'gff': str(p/'g.gff')}]
            result = list(builder.gene_generator(rows, [{'species':'s','chrom':'c'}], {'s':0}, flank=0))
            self.assertEqual(result[0]['input_ids'], [4,5,6])


if __name__ == '__main__':
    unittest.main()
