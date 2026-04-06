"""
    bilinear_interp(field, lat_grid, lon_grid, lat, lon)

Bilinear interpolation on a regular 2D grid. Assumes:
- `lat_grid` is a sorted vector of latitudes (dimension 1 of the field after setup)
- `lon_grid` is a sorted vector of longitudes (dimension 2 of the field after setup)
- `field` is size (length(lat_grid), length(lon_grid))
- Caller handles cyclic longitude wrapping before calling.

Returns interpolated value at (lat, lon).
"""
function bilinear_interp(field::AbstractMatrix, lat_grid::AbstractVector,
                         lon_grid::AbstractVector, lat::Real, lon::Real)
    # Find bracketing indices using binary search
    i = clamp(searchsortedfirst(lat_grid, lat) - 1, 1, length(lat_grid) - 1)
    j = clamp(searchsortedfirst(lon_grid, lon) - 1, 1, length(lon_grid) - 1)

    # Fractional positions
    t_lat = (lat - lat_grid[i]) / (lat_grid[i+1] - lat_grid[i])
    t_lon = (lon - lon_grid[j]) / (lon_grid[j+1] - lon_grid[j])

    # Clamp to [0,1] for edge cases
    t_lat = clamp(t_lat, 0.0, 1.0)
    t_lon = clamp(t_lon, 0.0, 1.0)

    # Bilinear blend
    return (1 - t_lat) * (1 - t_lon) * field[i, j] +
           (1 - t_lat) *      t_lon  * field[i, j+1] +
                t_lat  * (1 - t_lon) * field[i+1, j] +
                t_lat  *      t_lon  * field[i+1, j+1]
end

"""
    RegularGrid(lat_grid, lon_grid)

Precomputed structure for O(1) index lookup on a regular 2D grid.
"""
struct RegularGrid
    lat0::Float64
    lon0::Float64
    inv_dlat::Float64
    inv_dlon::Float64
    nlat::Int
    nlon::Int
end

function RegularGrid(lat_grid::AbstractVector, lon_grid::AbstractVector)
    dlat = lat_grid[2] - lat_grid[1]
    dlon = lon_grid[2] - lon_grid[1]
    RegularGrid(lat_grid[1], lon_grid[1], 1.0 / dlat, 1.0 / dlon,
                length(lat_grid), length(lon_grid))
end

"""
    bilinear_interp(field, grid::RegularGrid, lat, lon)

O(1) bilinear interpolation using precomputed regular grid spacing.
"""
function bilinear_interp(field::AbstractMatrix, grid::RegularGrid, lat::Real, lon::Real)
    # O(1) index computation via floor
    fi = (lat - grid.lat0) * grid.inv_dlat
    fj = (lon - grid.lon0) * grid.inv_dlon

    i = clamp(unsafe_trunc(Int, fi) + 1, 1, grid.nlat - 1)
    j = clamp(unsafe_trunc(Int, fj) + 1, 1, grid.nlon - 1)

    t_lat = clamp(fi - (i - 1), 0.0, 1.0)
    t_lon = clamp(fj - (j - 1), 0.0, 1.0)

    @inbounds begin
        return (1 - t_lat) * (1 - t_lon) * field[i, j] +
               (1 - t_lat) *      t_lon  * field[i, j+1] +
                    t_lat  * (1 - t_lon) * field[i+1, j] +
                    t_lat  *      t_lon  * field[i+1, j+1]
    end
end

"""
    RegularLookup(x, y)

Precomputed structure for O(1) 1D interpolation on regularly-spaced data.
"""
struct RegularLookup
    x0::Float64
    inv_dx::Float64
    n::Int
    y::Vector{Float64}
end

function RegularLookup(x::AbstractVector, y::AbstractVector)
    dx = x[2] - x[1]
    RegularLookup(x[1], 1.0 / dx, length(x), collect(Float64, y))
end

"""
    interp1_extrap(lut::RegularLookup, xq)

O(1) 1D linear interpolation with extrapolation using precomputed regular lookup.
"""
function interp1_extrap(lut::RegularLookup, xq::Real)
    fi = (xq - lut.x0) * lut.inv_dx
    lo = clamp(unsafe_trunc(Int, fi) + 1, 1, lut.n - 1)
    t = fi - (lo - 1)
    @inbounds return lut.y[lo] + t * (lut.y[lo+1] - lut.y[lo])
end

"""
    interp1_extrap(lut::RegularLookup, xq_vec::AbstractVector)

Vectorized O(1) 1D interpolation.
"""
function interp1_extrap(lut::RegularLookup, xq_vec::AbstractVector)
    return [interp1_extrap(lut, xq) for xq in xq_vec]
end

"""
    interp1_extrap(x, y, xq)

1D linear interpolation with extrapolation.
`x` must be sorted. Returns interpolated value at query point `xq`.
Extrapolates linearly beyond the data range.
"""
function interp1_extrap(x::AbstractVector, y::AbstractVector, xq::Real)
    n = length(x)

    if xq <= x[1]
        # Extrapolate below
        if n == 1
            return y[1]
        end
        t = (xq - x[1]) / (x[2] - x[1])
        return y[1] + t * (y[2] - y[1])
    elseif xq >= x[n]
        # Extrapolate above
        if n == 1
            return y[1]
        end
        t = (xq - x[n-1]) / (x[n] - x[n-1])
        return y[n-1] + t * (y[n] - y[n-1])
    else
        # Interpolate
        lo = searchsortedlast(x, xq)
        lo = clamp(lo, 1, n - 1)
        t = (xq - x[lo]) / (x[lo+1] - x[lo])
        return y[lo] + t * (y[lo+1] - y[lo])
    end
end

"""
    interp1_extrap(x, y, xq_vec::AbstractVector)

Vectorized version: interpolate at each point in xq_vec.
"""
function interp1_extrap(x::AbstractVector, y::AbstractVector, xq_vec::AbstractVector)
    return [interp1_extrap(x, y, xq) for xq in xq_vec]
end
