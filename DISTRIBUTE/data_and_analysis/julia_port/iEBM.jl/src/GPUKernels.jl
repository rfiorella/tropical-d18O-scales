"""
    GPUKernels.jl — GPU acceleration stubs (Phase 3).

Currently returns CPU backend for all operations.
"""

struct CPUBackend end

"""
    select_backend(config) -> CPUBackend

Select computation backend. Returns CPU for Phase 1.
"""
select_backend(config::RunConfig) = CPUBackend()

"""
    to_device(data, backend::CPUBackend)

Identity function for CPU backend.
"""
to_device(data, ::CPUBackend) = data
