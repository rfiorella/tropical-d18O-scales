"""
    TauBar.jl — Tau-bar integration (no decomposition).

Ported from Isotopes.taubar_noDecomposition() in attenuationMod_fxns.py (lines 184-556),
plus dist_to_tau_coords(), terrestrial_E_frac(), build_streamlines_df(),
and linear_idx_to_degrees() from accessory_fxns.py.
"""

using DataFrames

# --- Utility functions ---

"""
    linear_idx_to_degrees(idx_input, deg_grid) -> Vector

Convert 1-D indices [0, n] to degree values using y = mx + b.
Port of linear_idx_to_degrees from accessory_fxns.py.
"""
function linear_idx_to_degrees(idx_input::AbstractVector, deg_grid::AbstractVector)
    degmin = minimum(deg_grid)
    degmax = maximum(deg_grid)
    idxmin = 0
    idxmax = length(deg_grid)
    m_slope = (degmax - degmin) / (idxmax - idxmin)
    b_int = degmax - m_slope * idxmax
    return m_slope .* idx_input .+ b_int
end

"""
    dist_to_tau_coords(tau, nanindex, mu, wp, E0) -> (mu_taugrid, wp_taugrid, E0_taugrid, tausteps)

Resample from spatial coordinates to uniform tau coordinates.
Port of dist_to_tau_coords from accessory_fxns.py.
"""
function dist_to_tau_coords(tau::AbstractVector{T}, nanindex::BitVector,
                            mu::AbstractVector{T}, wp::AbstractVector{T},
                            E0::AbstractVector{T}) where {T}
    tau_valid = tau[nanindex]
    E0_valid = E0[nanindex]
    tau_max = maximum(tau_valid)
    tau_min = minimum(tau_valid)
    n_tausteps = length(tau_valid)

    tausteps = range(tau_min, tau_max, length=n_tausteps)

    # Linear interpolation from tau_valid onto uniform tausteps
    mu_taugrid = _interp1(collect(tausteps), tau_valid, mu)
    wp_taugrid = _interp1(collect(tausteps), tau_valid, wp)
    E0_taugrid = _interp1(collect(tausteps), tau_valid, E0_valid)

    return mu_taugrid, wp_taugrid, E0_taugrid, collect(tausteps)
end

"""
    _interp1(xq, x, y) -> Vector

1-D linear interpolation of (x, y) evaluated at query points xq.
Matches np.interp behavior (clamps outside range).
"""
function _interp1(xq::AbstractVector{T}, x::AbstractVector{T},
                  y::AbstractVector{T}) where {T}
    n = length(x)
    nq = length(xq)
    out = Vector{T}(undef, nq)

    for i in 1:nq
        xi = xq[i]
        if xi <= x[1]
            out[i] = y[1]
        elseif xi >= x[n]
            out[i] = y[n]
        else
            # Binary search for bracket
            lo = searchsortedlast(x, xi)
            lo = clamp(lo, 1, n - 1)
            t = (xi - x[lo]) / (x[lo + 1] - x[lo])
            out[i] = y[lo] + t * (y[lo + 1] - y[lo])
        end
    end
    return out
end

"""
    terrestrial_efrac(land_check, wp, coast_step, tau, E0, nanindex, n_samples)

Compute the terrestrial evaporation fraction and E-weighted tau-bar.
Port of terrestrial_E_frac from accessory_fxns.py.

Returns `(tau_wtd_mean, land_Esource_frac, E_terrFrac)`.
"""
function terrestrial_efrac(land_check::Real, wp::AbstractVector{T},
                           coast_step::Int, tau::AbstractVector{T},
                           E0::AbstractVector{T}, nanindex::BitVector,
                           n_samples::Int) where {T}
    # NOTE: Python uses np.sum(wp[:coast_step]) here, but wp is normalized by
    # np.trapz. This means sum(wp) > 1.0 in general, so E_terrFrac can exceed 1.0.
    # This matches the original Python behavior. See land_Esource_frac below which
    # correctly uses trapz/trapz for a bounded [0,1] fraction.
    E_terrFrac = coast_step > 0 ? sum(@view(wp[1:min(coast_step, length(wp))])) : zero(T)

    if land_check > T(0.6)
        land_wp = coast_step > 0 ? wp[1:min(coast_step, length(wp))] : T[]
        nsteps_land = length(land_wp)

        if nsteps_land == 0
            return zero(T), zero(T), E_terrFrac
        end

        wp_total = _trapz(wp)
        land_Esource_frac = wp_total > 0 ? _trapz(land_wp) / wp_total : zero(T)

        ns = min(n_samples, nsteps_land)

        # Compute sampling indices (matching np.linspace(..., dtype=int))
        samples_idx = round.(Int, range(1, nsteps_land, length=ns))

        tau_nan = tau[nanindex]
        E0_nan = E0[nanindex]
        ntau = length(tau_nan)

        this_wp = Vector{T}(undef, ns)
        this_taubar = Vector{T}(undef, ns)

        for (idx, si) in enumerate(samples_idx)
            this_wp[idx] = land_wp[si]
            if si <= ntau
                seg_tau = @view tau_nan[si:end]
                seg_E0 = @view E0_nan[si:end]
                denom = _trapz_exp(seg_E0, seg_tau)
                if denom > 0
                    numer = _trapz_tau_exp(seg_tau, seg_E0, seg_tau)
                    this_taubar[idx] = numer / denom
                else
                    this_taubar[idx] = zero(T)
                end
            else
                this_taubar[idx] = zero(T)
            end
        end

        total_wp = sum(this_wp)
        tau_wtd_mean = total_wp > 0 ? sum(this_taubar .* this_wp) / total_wp : zero(T)

        return tau_wtd_mean, land_Esource_frac, E_terrFrac
    else
        return zero(T), zero(T), E_terrFrac
    end
end

# Helper: trapz(E0 * exp(-tau))
function _trapz_exp(E0::AbstractVector{T}, tau::AbstractVector{T}) where {T}
    n = length(E0)
    n <= 1 && return zero(T)
    s = zero(T)
    @inbounds for i in 1:n-1
        s += (E0[i] * exp(-tau[i]) + E0[i+1] * exp(-tau[i+1])) / 2
    end
    return s
end

# Helper: trapz(tau * E0 * exp(-tau))
function _trapz_tau_exp(tau::AbstractVector{T}, E0::AbstractVector{T},
                        tau_full::AbstractVector{T}) where {T}
    n = length(tau)
    n <= 1 && return zero(T)
    s = zero(T)
    @inbounds for i in 1:n-1
        s += (tau[i] * E0[i] * exp(-tau[i]) + tau[i+1] * E0[i+1] * exp(-tau[i+1])) / 2
    end
    return s
end

"""
    compute_taubar_at_point(mu_taugrid, E0_taugrid, tausteps) -> Float64

The core tau-bar formula:
    tau_bar = trapz(E0 * (1/mu) * tau * exp(-tau)) / trapz(E0 * (1/mu) * exp(-tau))
"""
function compute_taubar_at_point(mu_taugrid::AbstractVector{T},
                                 E0_taugrid::AbstractVector{T},
                                 tausteps::AbstractVector{T}) where {T}
    n = length(tausteps)
    n <= 1 && return zero(T)

    numer = zero(T)
    denom = zero(T)
    @inbounds for i in 1:n-1
        w1 = E0_taugrid[i] / mu_taugrid[i] * exp(-tausteps[i])
        w2 = E0_taugrid[i+1] / mu_taugrid[i+1] * exp(-tausteps[i+1])
        numer += (w1 * tausteps[i] + w2 * tausteps[i+1]) / 2
        denom += (w1 + w2) / 2
    end
    return denom > 0 ? numer / denom : zero(T)
end

# --- Streamline DataFrame builder ---

"""
    build_streamlines_df!(df_rows, lat_grid, lon_grid, x, y, idx_counter,
                          result, nanindex, streamline_type;
                          save_streamline, coarsener, ignore_vars)

Collect streamline data into pre-allocated row vectors for later DataFrame construction.
Port of build_streamlines_df from accessory_fxns.py.

Returns a DataFrame of streamline rows for this point.
"""
function build_streamlines_df(lat_grid::AbstractVector, lon_grid::AbstractVector,
                              x::Int, y::Int, idx_counter::Int,
                              result::StreamlineResult{T}, nanindex::BitVector,
                              streamline_type::String;
                              save_streamline::Bool=false,
                              coarsener::Int=1) where {T}
    if save_streamline
        # Extract valid-index data
        nvalid = count(nanindex)
        dist_x = result.dist[nanindex]
        lat_x = result.lat_save[nanindex]
        lon_x = result.lon_save[nanindex]
        E0_x = result.E0[nanindex]
        Fmag_x = result.Fmag0[nanindex]
        P0_x = result.P0[nanindex]
        PminE_x = result.PminE[nanindex]
        tau_x = result.tau[nanindex]
        mu_x = result.mu
        wp_x = result.wp
        lfrac_x = result.lfrac0[nanindex]
        nx = coarsener
    else
        # Placeholder data
        dist_x = T[-9999]
        lat_x = T[-9999]
        lon_x = T[-9999]
        E0_x = T[-9999]
        Fmag_x = T[-9999]
        P0_x = T[-9999]
        PminE_x = T[-9999]
        tau_x = T[-9999]
        mu_x = T[-9999]
        wp_x = T[-9999]
        lfrac_x = T[-9999]
        nx = 1
    end

    # Subsample
    indices = 1:nx:length(dist_x)

    df = DataFrame(
        dist_km        = dist_x[indices],
        lat_streamline = lat_x[indices],
        lon_streamline = lon_x[indices],
        ET             = E0_x[indices],
        Fmag           = Fmag_x[indices],
        PRECT          = P0_x[indices],
        PminE          = PminE_x[indices],
        tau            = tau_x[indices],
        mu             = mu_x[min.(indices, length(mu_x))],
        wp             = wp_x[min.(indices, length(wp_x))],
        lfrac          = lfrac_x[indices],
        lat_sink       = fill(lat_grid[y], length(indices)),
        lon_sink       = fill(lon_grid[x], length(indices)),
        streamline_type = fill(streamline_type, length(indices)),
        idx            = fill(idx_counter, length(indices))
    )

    return df
end

# --- Main tau-bar computation ---

"""
    compute_taubar!(ds, config; iso_coords=nothing, runtype="isotopes_noEBM")

Main entry point for tau-bar computation without decomposition.
Port of Isotopes.taubar_noDecomposition() from attenuationMod_fxns.py.

Returns `(output_arrays::Dict, df_streamlines::DataFrame)`.
"""
function compute_taubar!(ds, config::RunConfig;
                         iso_coords::Union{Nothing,DataFrame}=nothing,
                         runtype::String="isotopes_noEBM")
    # Skip if not solving isotopes
    if occursin("noIsotopes", runtype)
        return Dict{String,Any}(), DataFrame()
    end

    iso = config.isotope
    compute_method = iso.compute_method
    streamline_collect = iso.collect_streamline_data
    land_only = iso.compute_land_only
    n_samples = iso.n_samples_taubar_land

    dx = iso.streamline_dx_km
    Dx = dx / 111.0
    taumax = iso.streamline_max_tau
    dmax = iso.streamline_max_dist_km
    Nmax = ceil(Int, dmax / dx)

    # Extract grid
    lat = Float64.(ds["lat"][:])
    lon = Float64.(ds["lon"][:])
    nlat = length(lat)
    nlon = length(lon)

    # Initialize arrays
    P, E, UQ, VQ, Fmag, lfrac = tau_array_initialize(ds, lat, lon, config)

    # Build interpolators
    Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = build_interpolators(
        lat, lon, E, P, UQ, VQ, lfrac)

    # Determine computation grid
    if compute_method == "bbox"
        lat_range = iso.bbox_lat_range
        lon_range = iso.bbox_lon_range
        res = iso.bbox_resolution

        # Build computation grid at specified resolution
        comp_lat = collect(Float64, range(-90, 90, step=res[1]))
        comp_lon = collect(Float64, range(0, 360, step=res[2]))

        lat_mask = (comp_lat .>= lat_range[1]) .& (comp_lat .<= lat_range[2])
        lon_mask = (comp_lon .>= lon_range[1]) .& (comp_lon .<= lon_range[2])

        A_slice = findall(lat_mask)
        B_slice = findall(lon_mask)
        grid_lat = comp_lat
        grid_lon = comp_lon
    else
        # coord_list mode uses dataset grid
        grid_lat = lat
        grid_lon = lon
        A_slice = Int[]
        B_slice = Int[]
    end

    # Output arrays
    tau_bar = fill(NaN, length(grid_lat), length(grid_lon))
    tau_bar_Ewtd_land = fill(NaN, length(grid_lat), length(grid_lon))
    moisture_dist_inland = fill(NaN, length(grid_lat), length(grid_lon))
    land_frac_of_streamline = fill(NaN, length(grid_lat), length(grid_lon))
    land_frac_of_evapSource = fill(NaN, length(grid_lat), length(grid_lon))

    # Pre-allocate streamline workspace
    result = StreamlineResult(Nmax)
    all_streamlines = DataFrame[]

    # Determine which points to save streamlines for
    xsave = Int[]
    ysave = Int[]
    if streamline_collect == "some" && iso_coords !== nothing
        _setup_streamline_saves!(xsave, ysave, iso_coords, grid_lat, grid_lon)
    end

    if compute_method == "bbox"
        n_cells = length(A_slice) * length(B_slice)
        idx_counter = 1

        for y_idx in A_slice
            for x_idx in B_slice
                lat0 = grid_lat[y_idx]
                lon0 = grid_lon[x_idx]

                # Check whether to save streamline
                save_streamline = _should_save_streamline(
                    streamline_collect, x_idx, y_idx, xsave, ysave)

                # Check if on land
                land_check = _nearest_landfrac(lfrac, lat, lon, lat0, lon0)
                if land_only && land_check == 0.0
                    moisture_dist_inland[y_idx, x_idx] = 0.0
                    idx_counter += 1
                    continue
                elseif land_check < 0.6
                    moisture_dist_inland[y_idx, x_idx] = 0.0
                end

                # March streamline
                march_streamline!(result, lat0, lon0, Nmax, taumax, Dx, dx,
                                  Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)

                # Terrestrial E fraction
                tau_wtd_mean, land_Esource_frac, E_terrFrac = terrestrial_efrac(
                    land_check, result.wp, result.coast_step,
                    result.tau, result.E0, result.nanindex, n_samples)

                tau_bar_Ewtd_land[y_idx, x_idx] = E_terrFrac * tau_wtd_mean
                moisture_dist_inland[y_idx, x_idx] = result.dist_inland
                max_dist = maximum(@view(result.dist[1:result.nsteps]))
                land_frac_of_streamline[y_idx, x_idx] = max_dist > 0 ? result.dist_inland / max_dist : 0.0
                land_frac_of_evapSource[y_idx, x_idx] = land_Esource_frac

                # Build streamline DataFrame
                df_sl = build_streamlines_df(
                    grid_lat, grid_lon, x_idx, y_idx, idx_counter,
                    result, result.nanindex, "full_model";
                    save_streamline=save_streamline,
                    coarsener=iso.streamline_save_coarsener)
                push!(all_streamlines, df_sl)

                # Resample to tau coordinates and compute tau-bar
                mu_taugrid, wp_taugrid, E0_taugrid, tausteps = dist_to_tau_coords(
                    result.tau, result.nanindex, result.mu, result.wp, result.E0)
                tau_bar[y_idx, x_idx] = compute_taubar_at_point(mu_taugrid, E0_taugrid, tausteps)

                # Progress
                pct = round(idx_counter / n_cells * 100, digits=2)
                println("  tau-bar: $pct% complete")
                idx_counter += 1
            end
        end

    elseif compute_method == "coord_list" && iso_coords !== nothing
        n_cells = nrow(iso_coords)
        idx_counter = 1

        for coord in 1:n_cells
            lat0 = Float64(iso_coords.lat[coord])
            lon0 = Float64(iso_coords.lon[coord])

            # Convert to 0-360
            if lon0 < 0
                lon0 = mod(lon0, 360.0)
            end

            save_streamline = _should_save_streamline_coord(
                streamline_collect, iso_coords, coord)

            # Find nearest grid indices
            y_idx = argmin(abs.(lat .- lat0))
            x_idx = argmin(abs.(lon .- lon0))

            land_check = lfrac[y_idx, x_idx]
            if land_only && land_check == 0.0
                moisture_dist_inland[y_idx, x_idx] = 0.0
                idx_counter += 1
                continue
            elseif land_check < 0.6
                moisture_dist_inland[y_idx, x_idx] = 0.0
            end

            march_streamline!(result, lat0, lon0, Nmax, taumax, Dx, dx,
                              Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)

            tau_wtd_mean, land_Esource_frac, E_terrFrac = terrestrial_efrac(
                land_check, result.wp, result.coast_step,
                result.tau, result.E0, result.nanindex, n_samples)

            tau_bar_Ewtd_land[y_idx, x_idx] = E_terrFrac * tau_wtd_mean
            moisture_dist_inland[y_idx, x_idx] = result.dist_inland
            max_dist = maximum(@view(result.dist[1:result.nsteps]))
            land_frac_of_streamline[y_idx, x_idx] = max_dist > 0 ? result.dist_inland / max_dist : 0.0
            land_frac_of_evapSource[y_idx, x_idx] = land_Esource_frac

            df_sl = build_streamlines_df(
                lat, lon, x_idx, y_idx, idx_counter,
                result, result.nanindex, "full_model";
                save_streamline=save_streamline,
                coarsener=iso.streamline_save_coarsener)
            push!(all_streamlines, df_sl)

            mu_taugrid, wp_taugrid, E0_taugrid, tausteps = dist_to_tau_coords(
                result.tau, result.nanindex, result.mu, result.wp, result.E0)
            tau_bar[y_idx, x_idx] = compute_taubar_at_point(mu_taugrid, E0_taugrid, tausteps)

            pct = round(idx_counter / n_cells * 100, digits=2)
            println("  tau-bar: $pct% complete")
            idx_counter += 1
        end
    end

    df_streamlines = isempty(all_streamlines) ? DataFrame() : vcat(all_streamlines...)

    output = Dict{String,Any}(
        "lat" => grid_lat,
        "lon" => grid_lon,
        "tau_bar" => tau_bar,
        "tau_bar_wtdLandEvap" => tau_bar_Ewtd_land,
        "moisture_dist_inland" => moisture_dist_inland,
        "streamline_frac_land" => land_frac_of_streamline,
        "Esource_frac_land" => land_frac_of_evapSource,
    )

    return output, df_streamlines
end

# --- Internal helpers ---

function _nearest_landfrac(lfrac::Matrix{Float64}, lat::Vector{Float64},
                           lon::Vector{Float64}, lat0::Real, lon0::Real)
    y = argmin(abs.(lat .- lat0))
    x = argmin(abs.(lon .- lon0))
    return lfrac[y, x]
end

function _setup_streamline_saves!(xsave::Vector{Int}, ysave::Vector{Int},
                                  iso_coords::DataFrame,
                                  grid_lat::AbstractVector, grid_lon::AbstractVector)
    if "collect_streamline" in names(iso_coords)
        collect_rows = iso_coords[iso_coords.collect_streamline .== "Y", :]
        for row in eachrow(collect_rows)
            lat_val = Float64(row.lat)
            lon_val = Float64(row.lon)
            if lon_val < 0
                lon_val = mod(lon_val, 360.0)
            end
            yi = argmin(abs.(grid_lat .- lat_val))
            xi = argmin(abs.(grid_lon .- lon_val))
            push!(ysave, yi)
            push!(xsave, xi)
        end
    end
end

function _should_save_streamline(collect_mode::String, x::Int, y::Int,
                                 xsave::Vector{Int}, ysave::Vector{Int})
    collect_mode == "all" && return true
    collect_mode == "none" && return false
    # "some": check if this point is in the save list
    for i in eachindex(xsave)
        if xsave[i] == x && ysave[i] == y
            return true
        end
    end
    return false
end

function _should_save_streamline_coord(collect_mode::String, iso_coords::DataFrame,
                                       coord_idx::Int)
    collect_mode == "all" && return true
    collect_mode == "none" && return false
    if "collect_streamline" in names(iso_coords)
        return iso_coords.collect_streamline[coord_idx] == "Y"
    end
    return false
end
