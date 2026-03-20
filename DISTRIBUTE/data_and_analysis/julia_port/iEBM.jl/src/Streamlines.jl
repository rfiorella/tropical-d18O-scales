"""
    Streamlines.jl — Core streamline marching kernel.

Ported from tau_streamline_point() in accessory_fxns.py (lines 56-177).
This is the performance-critical hot path.
"""

"""
Pre-allocated workspace for a single streamline march.
Reuse across grid cells to avoid allocations in the hot loop.
"""
mutable struct StreamlineResult{T<:AbstractFloat}
    tau::Vector{T}
    E0::Vector{T}
    P0::Vector{T}
    Fmag0::Vector{T}
    vq0::Vector{T}
    uq0::Vector{T}
    lfrac0::Vector{T}
    dist::Vector{T}
    lat_save::Vector{T}
    lon_save::Vector{T}
    PminE::Vector{T}
    # Derived quantities (computed post-march)
    mu::Vector{T}       # P/Fmag on valid indices
    wp::Vector{T}       # E*exp(-tau) / integral(E*exp(-tau)) on valid indices
    # Scalars
    coast_step::Int
    dist_inland::T
    nsteps::Int          # number of valid steps taken
    nanindex::BitVector  # mask for valid (non-zero tau or first point)
end

"""
    StreamlineResult(Nmax::Int, T=Float64)

Allocate a streamline workspace for at most `Nmax` steps.
"""
function StreamlineResult(Nmax::Int, ::Type{T}=Float64) where {T}
    return StreamlineResult{T}(
        zeros(T, Nmax), zeros(T, Nmax), zeros(T, Nmax),
        zeros(T, Nmax), zeros(T, Nmax), zeros(T, Nmax),
        zeros(T, Nmax), zeros(T, Nmax), zeros(T, Nmax),
        zeros(T, Nmax), zeros(T, Nmax),
        T[], T[],     # mu, wp — will be resized after march
        0, zero(T), 0,
        falses(Nmax)
    )
end

"""
    reset!(r::StreamlineResult)

Zero out all arrays and scalars for reuse.
"""
function reset!(r::StreamlineResult{T}) where {T}
    fill!(r.tau, zero(T))
    fill!(r.E0, zero(T))
    fill!(r.P0, zero(T))
    fill!(r.Fmag0, zero(T))
    fill!(r.vq0, zero(T))
    fill!(r.uq0, zero(T))
    fill!(r.lfrac0, zero(T))
    fill!(r.dist, zero(T))
    fill!(r.lat_save, zero(T))
    fill!(r.lon_save, zero(T))
    fill!(r.PminE, zero(T))
    fill!(r.nanindex, false)
    r.coast_step = 0
    r.dist_inland = zero(T)
    r.nsteps = 0
    return nothing
end

"""
    march_streamline!(result, lat0, lon0, Nmax, taumax, Dx, dx,
                      Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)

March upstream along moisture transport vectors from sink point (lat0, lon0).
Fills `result` in-place with tau profile, climatology along streamline, etc.

Port of `tau_streamline_point()` from accessory_fxns.py.

- `Dx` = dx / 111 (step size in degrees)
- `dx` = step size in km
"""
function march_streamline!(result::StreamlineResult{T},
                           lat0::Real, lon0::Real,
                           Nmax::Int, taumax::Real, Dx::Real, dx::Real,
                           Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit) where {T}
    reset!(result)

    lat = T(lat0)
    lon = T(lon0)

    # Interpolate initial point
    result.E0[1]     = Efit(lat, lon)
    result.P0[1]     = Pfit(lat, lon)
    result.vq0[1]    = vqfit(lat, lon)
    result.uq0[1]    = uqfit(lat, lon)
    result.Fmag0[1]  = Fmagfit(lat, lon)
    result.lfrac0[1] = lfracfit(lat, lon)
    result.PminE[1]  = result.P0[1] - result.E0[1]
    result.dist[1]   = zero(T)
    result.lat_save[1] = lat
    result.lon_save[1] = lon

    dist_inland = zero(T)
    coast_step = 0
    this_landmass = true
    step = 1  # 1-based: current valid position

    while result.tau[step] < taumax && step < Nmax
        # Direction from transport vectors
        dtheta = -result.vq0[step] / result.Fmag0[step] * T(Dx)

        # Pole clamping
        if lat + dtheta > T(89.5)
            dtheta = T(89.5) - lat
        end
        if lat + dtheta < T(-89.5)
            dtheta = T(-89.5) - lat
        end

        dphi = -result.uq0[step] / result.Fmag0[step] * T(Dx) / cosd(lat + dtheta / 2)

        # Distance in this step
        xdist = T(dx)
        ydist = T(dx) / cosd(lat + dtheta / 2)
        result.dist[step + 1] = result.dist[step] + sqrt(xdist^2 + ydist^2)

        lat1 = lat + dtheta
        lon1 = lon + dphi

        # Latitude wrapping past poles
        if lat1 < T(-90)
            lat1 = T(180) + lat1
            lon1 = lon1 - T(180)
        elseif lat1 > T(90)
            lat1 = T(180) - lat1
            lon1 = lon1 - T(180)
        end

        # Longitude wrapping to [0, 360]
        if lon1 < zero(T)
            lon1 = lon1 + T(360)
        elseif lon1 > T(360)
            lon1 = lon1 - T(360)
        end

        result.lat_save[step + 1] = lat1
        result.lon_save[step + 1] = lon1

        # Interpolate fields at new position
        result.vq0[step + 1]    = vqfit(lat1, lon1)
        result.uq0[step + 1]    = uqfit(lat1, lon1)
        result.E0[step + 1]     = Efit(lat1, lon1)
        result.Fmag0[step + 1]  = Fmagfit(lat1, lon1)
        result.P0[step + 1]     = Pfit(lat1, lon1)
        result.lfrac0[step + 1] = lfracfit(lat1, lon1)
        result.PminE[step + 1]  = result.P0[step + 1] - result.E0[step + 1]

        # Track distance over the same landmass
        if result.lfrac0[step + 1] > T(0.6) && this_landmass
            dist_inland = result.dist[step + 1]
        end
        if result.lfrac0[step + 1] <= T(0.6)
            this_landmass = false
            if coast_step == 0
                # Use step-1 so that coast_step indexes into the nanindex-compressed
                # wp array with the same semantics as Python's 0-based coast_step:
                # Python wp[:coast_step] gets coast_step elements;
                # Julia  wp[1:coast_step] also gets coast_step elements.
                coast_step = step - 1
            end
        end

        # Trapezoidal tau integration: dx*1000 converts km -> m
        result.tau[step + 1] = result.tau[step] +
            (dx * T(1000)) * (result.P0[step] / result.Fmag0[step] +
                              result.P0[step + 1] / result.Fmag0[step + 1]) / 2

        lat = lat1
        lon = lon1
        step += 1
    end

    # Handle case where air mass never crosses ocean
    # Use step-1 to match Python's 0-based coast_step semantics (see note above)
    if dist_inland > 0 && coast_step == 0
        coast_step = step - 1
    end

    result.coast_step = coast_step
    result.dist_inland = dist_inland
    result.nsteps = step

    # Build valid-index mask: tau > 0 or first point
    for k in 1:step
        result.nanindex[k] = (result.tau[k] > 0) || (k == 1)
    end

    # Interpolate E0 NaNs (rare but can occur)
    _interp_nans!(@view(result.E0[1:step]))

    # Compute derived quantities on valid indices
    valid_count = count(@view(result.nanindex[1:step]))
    result.mu = Vector{T}(undef, valid_count)
    result.wp = Vector{T}(undef, valid_count)

    idx = 0
    wp_sum = zero(T)
    for k in 1:step
        if result.nanindex[k]
            idx += 1
            result.mu[idx] = result.P0[k] / result.Fmag0[k]
            result.wp[idx] = result.E0[k] * exp(-result.tau[k])
            wp_sum += result.wp[idx]
        end
    end
    # Normalize wp (approximate integral normalization matching np.trapz behavior)
    if wp_sum > 0
        wp_integral = _trapz(@view(result.wp[1:valid_count]))
        if wp_integral > 0
            result.wp ./= wp_integral
        end
    end

    return nothing
end

"""
Simple trapezoidal integration (unit spacing), matching `np.trapz(y)`.
"""
function _trapz(y::AbstractVector{T}) where {T}
    n = length(y)
    n <= 1 && return zero(T)
    s = zero(T)
    @inbounds for i in 1:n-1
        s += (y[i] + y[i+1]) / 2
    end
    return s
end

"""
Trapezoidal integration with explicit x values, matching `np.trapz(y, x)`.
"""
function _trapz(y::AbstractVector{T}, x::AbstractVector{T}) where {T}
    n = length(y)
    n <= 1 && return zero(T)
    s = zero(T)
    @inbounds for i in 1:n-1
        s += (y[i] + y[i+1]) / 2 * (x[i+1] - x[i])
    end
    return s
end

"""
In-place linear interpolation to fill NaN values in a vector.
"""
function _interp_nans!(v::AbstractVector{T}) where {T}
    n = length(v)
    nan_present = false
    for i in 1:n
        if isnan(v[i])
            nan_present = true
            break
        end
    end
    nan_present || return nothing

    # Collect non-NaN indices and values
    good_idx = Int[]
    good_val = T[]
    for i in 1:n
        if !isnan(v[i])
            push!(good_idx, i)
            push!(good_val, v[i])
        end
    end
    isempty(good_idx) && return nothing

    # Linear interpolation
    for i in 1:n
        isnan(v[i]) || continue
        if i <= good_idx[1]
            v[i] = good_val[1]
        elseif i >= good_idx[end]
            v[i] = good_val[end]
        else
            # Find bracketing good indices
            lo = searchsortedlast(good_idx, i)
            hi = lo + 1
            t = (i - good_idx[lo]) / (good_idx[hi] - good_idx[lo])
            v[i] = good_val[lo] + t * (good_val[hi] - good_val[lo])
        end
    end
    return nothing
end
