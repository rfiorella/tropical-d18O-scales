"""
    EBM.jl — Energy Balance Model stubs.

The EBM functions are currently commented out in the Python source.
These stubs define the interface contracts for future implementation.
"""

"""
    solve_ebm!(ds, config) -> ds

Stub: solve the energy balance model. Currently a no-op.
"""
function solve_ebm!(ds, config::RunConfig)
    return ds
end

"""
    compute_efe!(ds, config) -> ds

Stub: compute the Energy Flux Equator. Currently a no-op.
"""
function compute_efe!(ds, config::RunConfig)
    return ds
end

"""
    compute_efpm!(ds, config) -> ds

Stub: compute the Energy Flux Prime Meridian. Currently a no-op.
"""
function compute_efpm!(ds, config::RunConfig)
    return ds
end
