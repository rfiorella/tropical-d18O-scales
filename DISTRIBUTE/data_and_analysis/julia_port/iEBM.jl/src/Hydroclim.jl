"""
    Hydroclim.jl — Orographic partitioning of precipitation fields.

Ported from Hydroclim.orog_partition() in attenuationMod_fxns.py (lines 43-166).
"""

"""
    compute_slope(elevation, lat, lon) -> Matrix{Float64}

Compute terrain slope magnitude from elevation grid using centered differences,
matching `np.gradient(elevation, lat, lon)`.
"""
function compute_slope(elevation::AbstractMatrix{<:Real},
                       lat::AbstractVector{<:Real},
                       lon::AbstractVector{<:Real})
    nlat, nlon = size(elevation)
    slope_y = zeros(Float64, nlat, nlon)
    slope_x = zeros(Float64, nlat, nlon)

    # Meridional gradient (along lat axis)
    for j in 1:nlon
        for i in 1:nlat
            if i == 1
                slope_y[i, j] = (elevation[2, j] - elevation[1, j]) / (lat[2] - lat[1])
            elseif i == nlat
                slope_y[i, j] = (elevation[nlat, j] - elevation[nlat-1, j]) / (lat[nlat] - lat[nlat-1])
            else
                slope_y[i, j] = (elevation[i+1, j] - elevation[i-1, j]) / (lat[i+1] - lat[i-1])
            end
        end
    end

    # Zonal gradient (along lon axis)
    for j in 1:nlon
        for i in 1:nlat
            if j == 1
                slope_x[i, j] = (elevation[i, 2] - elevation[i, 1]) / (lon[2] - lon[1])
            elseif j == nlon
                slope_x[i, j] = (elevation[i, nlon] - elevation[i, nlon-1]) / (lon[nlon] - lon[nlon-1])
            else
                slope_x[i, j] = (elevation[i, j+1] - elevation[i, j-1]) / (lon[j+1] - lon[j-1])
            end
        end
    end

    return sqrt.(slope_y .^ 2 .+ slope_x .^ 2)
end

"""
    build_orog_mask(slope, elevation, config) -> BitMatrix

Create orographic mask: true where slope >= threshold AND elevation > threshold.
"""
function build_orog_mask(slope::AbstractMatrix, elevation::AbstractMatrix,
                         config::RunConfig)
    return (slope .>= config.topo.slope_threshold_m_per_cell) .&
           (elevation .> config.topo.elev_threshold_m)
end

"""
    partition_field!(field, mask, landfrac, use_ocean)

Partition a 2-D field into orographic ("fixed") and background ("movable") components.

Returns `(fixed, movable)` where:
- `fixed` = anomaly relative to zonal mean, only where mask is true
- `movable` = field - fixed
"""
function partition_field(field::AbstractMatrix{T}, mask::AbstractMatrix{Bool},
                         landfrac::AbstractMatrix{T};
                         use_ocean::Bool=true) where {T<:AbstractFloat}
    nlat, nlon = size(field)

    # Compute zonal mean
    zonal_mean = zeros(T, nlat)
    for i in 1:nlat
        if use_ocean
            zonal_mean[i] = sum(@view field[i, :]) / nlon
        else
            count = 0
            s = zero(T)
            for j in 1:nlon
                if landfrac[i, j] > 0
                    s += field[i, j]
                    count += 1
                end
            end
            zonal_mean[i] = count > 0 ? s / count : zero(T)
        end
    end

    # Anomaly from zonal mean
    anomaly = similar(field)
    for j in 1:nlon
        for i in 1:nlat
            anomaly[i, j] = field[i, j] - zonal_mean[i]
        end
    end

    # Fixed = anomaly where mask is true, else 0
    fixed = similar(field)
    for j in 1:nlon
        for i in 1:nlat
            fixed[i, j] = mask[i, j] ? anomaly[i, j] : zero(T)
        end
    end

    movable = field .- fixed

    return fixed, movable
end
