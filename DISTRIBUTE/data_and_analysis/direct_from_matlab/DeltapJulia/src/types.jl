using Adapt

"""
    DeltapFields{T, M, V}

GPU-friendly struct holding all interpolation fields and precomputed grid parameters.

- `T`: Scalar type (Float32 or Float64)
- `M`: Matrix type (Array or GPU array)
- `V`: Vector type (Array or GPU vector)

All 2D fields are stored in transposed (nlat, nlon) layout matching the interpolation convention.
"""
struct DeltapFields{T, M <: AbstractMatrix{T}, V <: AbstractVector{T}}
    P::M
    UQ::M
    VQ::M
    E::M
    Tcond::M
    delta_e::M
    Plim::M
    lat0::T
    lon0::T
    inv_dlat::T
    inv_dlon::T
    n_lat::Int32
    n_lon::Int32
    alpha_T::V
    alpha_val::V
    alpha_Tmin::T
    alpha_inv_dT::T
end

Adapt.@adapt_structure DeltapFields
