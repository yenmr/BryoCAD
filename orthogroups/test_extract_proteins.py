"""Biological edge cases for strand, phase, isoform selection and QC."""
import csv
import importlib.util
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
spec=importlib.util.spec_from_file_location('extract',Path(__file__).with_name('extract_proteins_from_gff.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class ProteinTests(unittest.TestCase):
    def tx(self,parts,strand='+'):
        return {'seqid':'chr','strand':strand,'features':{'cds':[('chr',a,b,strand,str(p)) for a,b,p in parts]}}
    def test_split_codon_plus(self):
        aa,*_=m.translate_cds(self.tx([(0,4,0),(7,12,2)]),{'chr':'ATGACCCAATAA'})
        self.assertEqual(aa,'MK')
    def test_split_codon_minus(self):
        seq='ATGACCCAATAA'.translate(m.COMPLEMENT)[::-1]
        aa,*_=m.translate_cds(self.tx([(8,12,0),(0,5,2)],'-'),{'chr':seq})
        self.assertEqual(aa,'MK')
    def test_partial_initial_phase(self):
        aa,warnings,phase,rem=m.translate_cds(self.tx([(0,11,1)]),{'chr':'AATGAAATAAC'})
        self.assertEqual(aa,'MK');self.assertEqual(phase,1);self.assertEqual(rem,1)
        self.assertEqual(warnings,['partial_5prime','partial_3prime'])
    def test_bad_phase_rejected(self):
        with self.assertRaisesRegex(ValueError,'phase'):
            m.translate_cds(self.tx([(0,4,0),(7,12,0)]),{'chr':'ATGACCCAATAA'})
    def test_internal_stop_rejected(self):
        with self.assertRaisesRegex(ValueError,'internal_stop'):
            m.translate_cds(self.tx([(0,12,0)]),{'chr':'ATGTAAAAATAA'})
    def test_ambiguous_codon(self):
        aa,*_=m.translate_cds(self.tx([(0,9,0)]),{'chr':'ATGNNNTAA'})
        self.assertEqual(aa,'MX')
    def test_transcript_length_not_genomic_span(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=Path(tmp);(r/'g.fa').write_text('>chr\nATGAAATAA'+'A'*91+'\n')
            lines=[('gene',1,100,'.','ID=g'),('mRNA',1,100,'.','ID=wide;Parent=g'),
                ('exon',1,9,'.','Parent=wide'),('exon',99,100,'.','Parent=wide'),('CDS',1,9,'0','Parent=wide'),
                ('mRNA',1,20,'.','ID=long;Parent=g'),('exon',1,20,'.','Parent=long'),('CDS',1,9,'0','Parent=long')]
            (r/'g.gff').write_text('\n'.join(f'chr\t.\t{k}\t{a}\t{b}\t.\t+\t{phase}\t{attr}' for k,a,b,phase,attr in lines)+'\n')
            a=Namespace(fasta=r/'g.fa',gff=r/'g.gff',output=r/'p.fa',report=None,length_metric='transcript',id_suffix='|test',min_len=0,max_xfrac=.3,nterm_x=20)
            self.assertEqual(m.extract(a),1)
            with (r/'p.fa.selection.tsv').open() as f:row=next(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(row['transcript_id'],'long');self.assertEqual((r/'p.fa').read_text(),'>g|test\nMK\n')

if __name__=='__main__':unittest.main()
