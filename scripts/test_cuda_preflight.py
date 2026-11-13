from types import SimpleNamespace
import unittest
from cuda_preflight import check_cuda_runtime
class Value:
    def __add__(self,other):return self
    def sum(self):return self
    def item(self):return 36
class CudaTests(unittest.TestCase):
    def torch(self):
        return SimpleNamespace(version=SimpleNamespace(cuda='test'),float32='float32',arange=lambda *a,**k:Value(),cuda=SimpleNamespace(is_available=lambda:True,get_device_name=lambda i:'test GPU',get_device_capability=lambda i:(7,0),get_arch_list=lambda:['sm_70'],device_count=lambda:1,synchronize=lambda i:None))
    def test_visible_and_working(self):
        report,error=check_cuda_runtime(self.torch());self.assertIsNone(error);self.assertTrue(report['cuda_basic_operation_passed'])
    def test_visible_but_kernel_unavailable(self):
        t=self.torch()
        def fail(*a,**k):raise RuntimeError('no kernel image is available')
        t.arange=fail
        report,error=check_cuda_runtime(t);self.assertTrue(report['cuda_available']);self.assertFalse(report['cuda_basic_operation_passed']);self.assertIn('no kernel image',error)
    def test_no_gpu(self):
        t=self.torch();t.cuda.is_available=lambda:False
        report,error=check_cuda_runtime(t);self.assertFalse(report['cuda_available']);self.assertIsNotNone(error)
if __name__=='__main__':unittest.main()
