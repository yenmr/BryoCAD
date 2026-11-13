"""Distinguish visible CUDA devices from a working PyTorch CUDA build."""
def check_cuda_runtime(torch):
    result = {'cuda_available': torch.cuda.is_available(), 'cuda': torch.version.cuda}
    if not result['cuda_available']:
        return result, 'No CUDA GPU is visible'
    try:
        result['gpu'] = torch.cuda.get_device_name(0)
        result['compute_capability'] = list(torch.cuda.get_device_capability(0))
        result['compiled_architectures'] = torch.cuda.get_arch_list()
        result['visible_gpu_count'] = torch.cuda.device_count()
        # CUDA availability alone does not establish that any kernel can run.
        values = torch.arange(8, device='cuda:0', dtype=torch.float32)
        total = (values + 1).sum().item()
        torch.cuda.synchronize(0)
        if total != 36:
            raise RuntimeError(f'Unexpected CUDA probe result: {total}')
        result['cuda_basic_operation_passed'] = True
        return result, None
    except Exception as exc:
        result['cuda_basic_operation_passed'] = False
        return result, 'Basic CUDA operation failed with the selected PyTorch build: ' + str(exc)
