"""VCF serialization integration test; requires PyVCF3 (no model or GPU)."""
import math
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from test_snp_scoring import functions
try:
    import vcf
except ImportError:
    vcf = None

@unittest.skipIf(vcf is None, 'Install PyVCF3 to test actual VCF serialization')
class VcfIntegration(unittest.TestCase):
    def test_roundtrip_and_header_only(self):
        ns=functions({'zero_shot_score_vcf'}, {'vcf':vcf,'np':SimpleNamespace(log=math.log),
            'logging':SimpleNamespace(info=lambda *a:None)})
        with TemporaryDirectory() as d:
            source=Path(d)/'in.vcf'; output=Path(d)/'out.vcf'
            source.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t300\t.\tA\tC\t.\tPASS\t.\nchr1\t400\t.\tAT\tA\t.\tPASS\t.\n')
            args=SimpleNamespace(inputVCF=str(source),output=str(output))
            ns['zero_shot_score_vcf'](args,[0],[[.2,.4,.3,.1]])
            with output.open() as stream:
                reader=vcf.Reader(stream);records=list(reader)
            self.assertEqual(len(records),1)
            self.assertEqual(reader.infos['plantCAD_zero_shot'].num,-1)
            self.assertAlmostEqual(float(records[0].INFO['plantCAD_zero_shot'][0]),math.log(2))
            ns['zero_shot_score_vcf'](args,[],[])
            with output.open() as stream:
                self.assertEqual(list(vcf.Reader(stream)),[])

if __name__=='__main__':unittest.main()
