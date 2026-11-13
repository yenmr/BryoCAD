"""Dependency-free regression tests of VCF writing and context extraction."""
import ast
import math
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import unittest

SOURCE = Path(__file__).with_name('zero_shot_score.py').read_text()

def functions(names, namespace):
    tree = ast.parse(SOURCE)
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    exec(compile(ast.Module(body=nodes, type_ignores=[]), '<scorer functions>', 'exec'), namespace)
    return namespace

class ScorerTests(unittest.TestCase):
    def window(self):
        return functions({'centered_variant_window'}, {})['centered_variant_window']

    def test_short_contigs_keep_target_at_center(self):
        window = self.window()
        for length in [1, 10, 100, 511, 512, 513, 1024]:
            contig = ('ACGT' * (length // 4 + 1))[:length]
            for pos in range(length):
                result = window(contig, pos, 512, 255)
                self.assertEqual(len(result), 512)
                self.assertEqual(result[255], contig[pos])
                for offset, base in enumerate(result):
                    coordinate = pos - 255 + offset
                    self.assertEqual(base, contig[coordinate] if 0 <= coordinate < length else 'N')

    def test_long_contig_windows_identical_to_original(self):
        window = self.window(); contig = 'ACGT' * 400
        for pos in range(len(contig)):
            old = (contig[:pos + 257].rjust(512, 'N') if pos < 255
                   else contig[pos - 255:pos + 257].ljust(512, 'N'))
            self.assertEqual(window(contig, pos, 512, 255), old)

    def test_invalid_coordinates_rejected(self):
        for pos in [-1, 4]:
            with self.assertRaises(ValueError):self.window()('ACGT', pos, 512, 255)

    def write(self, records, indices, probs):
        output = []; headers = []
        class Reader:
            infos = {}
            def __iter__(self):return iter(records)
        class Writer:
            def __init__(self, handle, reader):
                handle.close(); headers.append(dict(reader.infos))
            def write_record(self, record):output.append(record)
            def close(self):pass
        ns = functions({'zero_shot_score_vcf'}, {
            'vcf': SimpleNamespace(Reader=lambda **kw:Reader(), Writer=Writer,
                                  parser=SimpleNamespace(_vcf_metadata_parser=lambda:SimpleNamespace(read_info=lambda line:('plantCAD_zero_shot',line)))),
            'logging': SimpleNamespace(info=lambda *a:None), 'np': SimpleNamespace(log=math.log)})
        with TemporaryDirectory() as tmp:
            ns['zero_shot_score_vcf'](SimpleNamespace(inputVCF='toy.vcf',output=str(Path(tmp)/'out.vcf')), indices, probs)
        self.assertIn('plantCAD_zero_shot',headers[0])
        return output

    @staticmethod
    def snp():
        return SimpleNamespace(REF='A', ALT=[SimpleNamespace(type='SNV',sequence='C')], INFO={})

    def test_trailing_non_snv_no_crash(self):
        result=self.write([self.snp(), SimpleNamespace(REF='A',ALT=[],INFO={})],[0],[[.2,.4,.3,.1]])
        self.assertEqual(len(result),1)
        self.assertAlmostEqual(float(result[0].INFO['plantCAD_zero_shot']),math.log(.4/.2))

    def test_empty_and_no_snv_inputs(self):
        self.assertEqual(self.write([],[],[]),[])
        self.assertEqual(self.write([SimpleNamespace(REF='A',ALT=[],INFO={})],[],[]),[])

    def test_multiple_scores_keep_formula_and_order(self):
        result=self.write([self.snp(),self.snp()],[0,1],[[.1,.2,.3,.4],[.4,.2,.3,.1]])
        self.assertEqual([x.INFO['plantCAD_zero_shot'] for x in result], [str(math.log(2)),str(math.log(.5))])

    def test_empty_input_skips_model_loading(self):
        called=[]
        ns=functions({'main'}, {'parse_args':lambda:SimpleNamespace(inputVCF='empty.vcf',output='out.vcf'),
            'seq_from_vcf':lambda a:([],[]), 'zero_shot_score_vcf':lambda *a:called.append(a),
            'load_model_and_tokenizer':lambda *a:self.fail('Empty input loaded model'),
            'logging':SimpleNamespace(basicConfig=lambda **k:None,INFO=20,info=lambda *a:None)})
        ns['main']();self.assertEqual(len(called),1)

if __name__ == '__main__':unittest.main()
