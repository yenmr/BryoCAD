"""Verify wrapper dispatch under the shallow-copy semantics used by DataParallel."""
import ast
from pathlib import Path
from types import SimpleNamespace
import unittest

class ReplicaTests(unittest.TestCase):
    def loader(self):
        tree=ast.parse(Path(__file__).with_name('lora_fine_tune.py').read_text())
        node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='load_base_model')
        class Model:
            def __init__(self):self.device='cuda:0'
            def forward(self,input_ids=None,labels=None):
                if input_ids != self.device:raise RuntimeError('Wrong replica dispatch')
                return self.device,labels
        factory=SimpleNamespace(from_pretrained=lambda *a,**k:Model())
        ns={'Optional':__import__('typing').Optional,'AutoModelForSequenceClassification':factory}
        exec(compile(ast.Module(body=[node],type_ignores=[]),'<load_base_model>','exec'),ns)
        return ns['load_base_model']
    def test_replica_uses_own_model_and_repeated_load_is_safe(self):
        load=self.loader()
        first=load('unused','classification',2)
        second=load('unused','classification',2)
        replica=object.__new__(type(first));replica.__dict__=first.__dict__.copy();replica.device='cuda:1'
        self.assertEqual(first(input_ids='cuda:0') if callable(first) else first.forward(input_ids='cuda:0'),('cuda:0',None))
        self.assertEqual(replica.forward(input_ids='cuda:1',attention_mask=None),('cuda:1',None))
        self.assertEqual(second.forward(input_ids='cuda:0',labels=1),('cuda:0',1))
    def test_multilabel_conversion_is_per_instance(self):
        load=self.loader();multi=load('unused','multi_label',2);single=load('unused','classification',2)
        label=SimpleNamespace(float=lambda:'converted')
        self.assertEqual(multi.forward(input_ids='cuda:0',labels=label),('cuda:0','converted'))
        self.assertIs(single.forward(input_ids='cuda:0',labels=label)[1],label)

if __name__=='__main__':unittest.main()
