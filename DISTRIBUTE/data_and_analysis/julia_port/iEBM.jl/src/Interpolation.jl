"""
    Interpolation.jl — Allocation-free bilinear interpolation on a regular grid.

Ported from tau_interpolatorInitialize() in accessory_fxns.py (lines 227-253).
Replaces scipy.interpolate.RegularGridInterpolator with method='linear'.
"""

"""
    RegularGridInterp{T}

Bilinear interpolator on a regular (lat, lon) grid with periodic longitude
wraparound. The longitude axis is padded by one cell on each side so that
interpolation near the 0/360 boundary works correctly.

Constructor pads the data; the callable `(interp)(lat, lon)` performs
allocation-free bilinear interpolation.
"""
struct RegularGridInterp{T<:AbstractFloat}
    lat::Vector{T}      # original lat grid (sorted ascending)
    lon_padded::Vector{T}  # lon grid with wraparound padding
    data::Matrix{T}      # [nlat, nlon_padded]

    # Precomputed inverse spacings for speed
    inv_dlat::T
    inv_dlon::T
end

"""
    RegularGridInterp(lat, lon, data)

Build a bilinear interpolator. `data` must be `[nlat, nlon]`.
Longitude is padded with one wraparound cell on each side, matching the Python:
    lon_test = [lon[-1]-360, lon..., lon[0]+360]
"""
function RegularGridInterp(lat::AbstractVector{T}, lon::AbstractVector{T},
                           data::AbstractMatrix{T}) where {T<:AbstractFloat}
    nlat, nlon = size(data)
    @assert length(lat) == nlat
    @assert length(lon) == nlon

    # Pad longitude with wraparound
    lon_padded = Vector{T}(undef, nlon + 2)
    lon_padded[1] = lon[end] - T(360)
    lon_padded[2:end-1] .= lon
    lon_padded[end] = lon[1] + T(360)

    # Pad data: wrap last column to front, first column to back
    data_padded = Matrix{T}(undef, nlat, nlon + 2)
    data_padded[:, 1]        .= @view data[:, end]
    data_padded[:, 2:end-1]  .= data
    data_padded[:, end]      .= @view data[:, 1]

    dlat = lat[2] - lat[1]
    dlon = lon[2] - lon[1]

    return RegularGridInterp{T}(
        Vector{T}(lat), lon_padded, data_padded,
        one(T) / dlat, one(T) / dlon
    )
end

"""
    (interp::RegularGridInterp)(lat, lon) -> value

Evaluate bilinear interpolation at (lat, lon). Allocation-free.
"""
function (interp::RegularGridInterp{T})(lat::Real, lon::Real) where {T}
    latv = T(lat)
    lonv = T(lon)

    lat_grid = interp.lat
    lon_grid = interp.lon_padded
    data = interp.data

    # Find lat index (1-based, clamped)
    fi = (latv - lat_grid[1]) * interp.inv_dlat + one(T)
    i = clamp(floor(Int, fi), 1, length(lat_grid) - 1)
    t_lat = clamp(fi - T(i), zero(T), one(T))

    # Find lon index (1-based)
    fj = (lonv - lon_grid[1]) * interp.inv_dlon + one(T)
    j = clamp(floor(Int, fj), 1, length(lon_grid) - 1)
    t_lon = clamp(fj - T(j), zero(T), one(T))

    # Bilinear interpolation
    v00 = data[i,     j]
    v10 = data[i + 1, j]
    v01 = data[i,     j + 1]
    v11 = data[i + 1, j + 1]

    return (one(T) - t_lat) * ((one(T) - t_lon) * v00 + t_lon * v01) +
           t_lat            * ((one(T) - t_lon) * v10 + t_lon * v11)
end

"""
    build_interpolators(lat, lon, E, P, UQ, VQ, lfrac)

Build all six interpolators needed for streamline marching. Matches
`tau_interpolatorInitialize` from accessory_fxns.py.

Returns `(Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit)`.
"""
function build_interpolators(lat::AbstractVector{T}, lon::AbstractVector{T},
                             E::AbstractMatrix{T}, P::AbstractMatrix{T},
                             UQ::AbstractMatrix{T}, VQ::AbstractMatrix{T},
                             lfrac::AbstractMatrix{T}) where {T<:AbstractFloat}
    Fmag = sqrt.(UQ .^ 2 .+ VQ .^ 2)

    Efit     = RegularGridInterp(lat, lon, E)
    Pfit     = RegularGridInterp(lat, lon, P)
    uqfit    = RegularGridInterp(lat, lon, UQ)
    vqfit    = RegularGridInterp(lat, lon, VQ)
    Fmagfit  = RegularGridInterp(lat, lon, Fmag)
    lfracfit = RegularGridInterp(lat, lon, lfrac)

    return Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit
end
