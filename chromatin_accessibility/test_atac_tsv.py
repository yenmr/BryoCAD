import tempfile,unittest
from pathlib import Path
from atac_tsv import read_atac_tsv
class TsvTests(unittest.TestCase):
    def read(self,text):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.tsv';p.write_text(text);return read_atac_tsv(p)
    def test_optional_labels_and_metadata(self):
        rows=self.read('sequence\tlabel\tchromosome\n'+'acgt'*128+'\t1\tchr7\n')
        self.assertEqual(rows[0]['sequence'],'ACGT'*128);self.assertEqual(rows[0]['label'],'1');self.assertEqual(rows[0]['chromosome'],'chr7')
        self.assertEqual(len(self.read('sequence\n'+'A'*512+'\n')),1)
    def test_wrong_length(self):
        with self.assertRaises(ValueError):self.read('sequence\nACGT\n')
    def test_bad_label(self):
        with self.assertRaises(ValueError):self.read('sequence\tlabel\n'+'A'*512+'\t2\n')
    def test_bad_base(self):
        with self.assertRaises(ValueError):self.read('sequence\n'+'A'*511+'X\n')
    def test_missing_fields(self):
        for text in ['sequence\n','label\n1\n','sequence\tlabel\n'+'A'*512+'\n']:
            with self.assertRaises(ValueError):self.read(text)
if __name__=='__main__':unittest.main()
