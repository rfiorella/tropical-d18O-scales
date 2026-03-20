"""
    Decomposition.jl — Attribution analysis (tau-bar decomposition).

Ported from Isotopes.taubar_withDecomposition() in attenuationMod_fxns.py (lines 558-1434),
plus tau_spatial_decompose(), dist_to_Dtaubar_frac_LocalEvap(),
dist_to_Dtaubar_frac_streamline(), local_plus_regional_tauStreamline(),
tau_decompose_ELPath() from accessory_fxns.py.
"""

using DataFrames

# --- Spatial decomposition ---

"""
    spatial_decompose(local_threshold, dist, tau, E0, nanindex,
                      disti, taui, E0i, nanindexi; local_eff_mode=false)

Decompose tau-bar change into local vs regional components.
Port of tau_spatial_decompose from accessory_fxns.py.

If `local_eff_mode=true`, returns just the local effect scalar.
Otherwise returns `(regional_effect, local_effect)`.
"""
function spatial_decompose(local_threshold::Real,
                           dist::AbstractVector{T}, tau::AbstractVector{T},
                           E0::AbstractVector{T}, nanindex::BitVector,
                           disti::AbstractVector{T}, taui::AbstractVector{T},
                           E0i::AbstractVector{T}, nanindexi::BitVector;
                           local_eff_mode::Bool=false) where {T}
    # Find local/remote indices
    dist_valid = dist[nanindex]
    disti_valid = disti[nanindexi]

    local_idx_old = findall(disti_valid .<= local_threshold)
    local_idx_new = findall(dist_valid .<= local_threshold)

    nandxmax = findlast(nanindex)
    nandxmaxi = findlast(nanindexi)

    # Old local vs regional
    max_local_old = isempty(local_idx_old) ? 0 : maximum(local_idx_old)
    max_local_new = isempty(local_idx_new) ? 0 : maximum(local_idx_new)

    tau_old_local = taui[local_idx_old]
    tau_old_regional = max_local_old < nandxmaxi ? taui[max_local_old+1:nandxmaxi] : T[]
    E0_old_local = E0i[local_idx_old]
    E0_old_regional = max_local_old < nandxmaxi ? E0i[max_local_old+1:nandxmaxi] : T[]

    tau_new_local = tau[local_idx_new]
    tau_new_regional = max_local_new < nandxmax ? tau[max_local_new+1:nandxmax] : T[]
    E0_new_local = E0[local_idx_new]
    E0_new_regional = max_local_new < nandxmax ? E0[max_local_new+1:nandxmax] : T[]

    # Composite: local constant = old local + new regional
    tau_localConst = vcat(tau_old_local, tau_new_regional)
    E0_localConst = vcat(E0_old_local, E0_new_regional)
    tau_regionalConst = vcat(tau_new_local, tau_old_regional)
    E0_regionalConst = vcat(E0_new_local, E0_old_regional)

    taubar_localConst = _taubar_simple(tau_localConst, E0_localConst)
    taubar_regionalConst = _taubar_simple(tau_regionalConst, E0_regionalConst)
    taubar_old = _taubar_simple(taui[nanindexi], E0i[nanindexi])

    if local_eff_mode
        return taubar_regionalConst - taubar_old
    else
        regional_effect = taubar_localConst - taubar_old
        local_effect = taubar_regionalConst - taubar_old
        return regional_effect, local_effect
    end
end

"""
Simple tau-bar = trapz(tau * E0 * exp(-tau)) / trapz(E0 * exp(-tau))
"""
function _taubar_simple(tau::AbstractVector{T}, E0::AbstractVector{T}) where {T}
    n = length(tau)
    n <= 1 && return zero(T)
    numer = zero(T)
    denom = zero(T)
    @inbounds for i in 1:n-1
        w1 = E0[i] * exp(-tau[i])
        w2 = E0[i+1] * exp(-tau[i+1])
        numer += (tau[i] * w1 + tau[i+1] * w2) / 2
        denom += (w1 + w2) / 2
    end
    return denom > 0 ? numer / denom : zero(T)
end

# --- Composite streamline construction ---

"""
    local_plus_regional_streamline(local_threshold, dx,
        P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
        P0, Fmag0, E0, dist, nanindex, lfrac0)

Construct composite streamlines for upwind effect decomposition.
Port of local_plus_regional_tauStreamline from accessory_fxns.py.

Returns (tau_localConst, E0_localConst, tau_regionalConst, E0_regionalConst,
         Fmag0_localConst, Fmag0_regionalConst, P0_localConst, P0_regionalConst,
         dist_localConst, dist_regionalConst).
"""
function local_plus_regional_streamline(local_threshold::Real, dx::Real,
        P0i::AbstractVector{T}, Fmag0i::AbstractVector{T},
        E0i::AbstractVector{T}, disti::AbstractVector{T},
        nanindexi::BitVector, lfrac0i::AbstractVector{T},
        P0::AbstractVector{T}, Fmag0::AbstractVector{T},
        E0::AbstractVector{T}, dist::AbstractVector{T},
        nanindex::BitVector, lfrac0::AbstractVector{T}) where {T}

    disti_valid = disti[nanindexi]
    dist_valid = dist[nanindex]

    local_idx_old = findall(disti_valid .<= local_threshold)
    local_idx_new = findall(dist_valid .<= local_threshold)
    nandxmax = findlast(nanindex)
    nandxmaxi = findlast(nanindexi)
    max_local_old = isempty(local_idx_old) ? 0 : maximum(local_idx_old)
    max_local_new = isempty(local_idx_new) ? 0 : maximum(local_idx_new)

    # Extract sub-arrays
    P0_old_local = P0i[local_idx_old]
    P0_old_regional = max_local_old < nandxmaxi ? P0i[max_local_old+1:nandxmaxi] : T[]
    Fmag0_old_local = Fmag0i[local_idx_old]
    Fmag0_old_regional = max_local_old < nandxmaxi ? Fmag0i[max_local_old+1:nandxmaxi] : T[]
    E0_old_local = E0i[local_idx_old]
    E0_old_regional = max_local_old < nandxmaxi ? E0i[max_local_old+1:nandxmaxi] : T[]
    dist_old_local = disti[local_idx_old]
    dist_old_regional = max_local_old < nandxmaxi ? disti[max_local_old+1:nandxmaxi] : T[]

    P0_new_local = P0[local_idx_new]
    P0_new_regional = max_local_new < nandxmax ? P0[max_local_new+1:nandxmax] : T[]
    Fmag0_new_local = Fmag0[local_idx_new]
    Fmag0_new_regional = max_local_new < nandxmax ? Fmag0[max_local_new+1:nandxmax] : T[]
    E0_new_local = E0[local_idx_new]
    E0_new_regional = max_local_new < nandxmax ? E0[max_local_new+1:nandxmax] : T[]
    dist_new_local = dist[local_idx_new]
    dist_new_regional = max_local_new < nandxmax ? dist[max_local_new+1:nandxmax] : T[]

    # Build composites
    P0_localConst = vcat(P0_old_local, P0_new_regional)
    P0_regionalConst = vcat(P0_new_local, P0_old_regional)
    Fmag0_localConst = vcat(Fmag0_old_local, Fmag0_new_regional)
    Fmag0_regionalConst = vcat(Fmag0_new_local, Fmag0_old_regional)
    E0_localConst = vcat(E0_old_local, E0_new_regional)
    E0_regionalConst = vcat(E0_new_local, E0_old_regional)
    dist_localConst = vcat(dist_old_local, dist_new_regional)
    dist_regionalConst = vcat(dist_new_local, dist_old_regional)

    # Integrate tau for each composite
    tau_localConst = _integrate_tau(P0_localConst, Fmag0_localConst, dx)
    tau_regionalConst = _integrate_tau(P0_regionalConst, Fmag0_regionalConst, dx)

    return (tau_localConst, E0_localConst, tau_regionalConst, E0_regionalConst,
            Fmag0_localConst, Fmag0_regionalConst, P0_localConst, P0_regionalConst,
            dist_localConst, dist_regionalConst)
end

"""Trapezoidal tau integration from P and Fmag arrays."""
function _integrate_tau(P0::AbstractVector{T}, Fmag0::AbstractVector{T}, dx::Real) where {T}
    n = length(P0)
    tau = zeros(T, n)
    @inbounds for i in 1:n-1
        tau[i+1] = tau[i] + (dx * T(1000)) * (P0[i] / Fmag0[i] + P0[i+1] / Fmag0[i+1]) / 2
    end
    return tau
end

# --- Binary search for fraction distance ---

"""
    find_fraction_distance_localevap(Dtau_bar_val, dist, tau, E0, nanindex,
        disti, taui, E0i, nanindexi; dtau_frac=0.75, perc_error=1.0)

Binary search for distance where `dtau_frac` of Dtau_bar is explained
by local evaporation effect.
Port of dist_to_Dtaubar_frac_LocalEvap from accessory_fxns.py.
"""
function find_fraction_distance_localevap(Dtau_bar_val::T,
        dist::AbstractVector{T}, tau::AbstractVector{T},
        E0::AbstractVector{T}, nanindex::BitVector,
        disti::AbstractVector{T}, taui::AbstractVector{T},
        E0i::AbstractVector{T}, nanindexi::BitVector;
        dtau_frac::Float64=0.75, perc_error::Float64=1.0) where {T}

    this_Dtau_frac = Dtau_bar_val * dtau_frac
    min_Dtau_frac = minimum([this_Dtau_frac * (1 - perc_error/100),
                             this_Dtau_frac * (1 + perc_error/100)])
    max_Dtau_frac = maximum([this_Dtau_frac * (1 - perc_error/100),
                             this_Dtau_frac * (1 + perc_error/100)])

    mindist = 0
    maxdist = Int(round(maximum(filter(!isnan, dist))))
    localEff = T(-1e3)
    testdistance = 0
    searchiter = 0
    binary_search_mode = true
    brute_force_step = 5
    searchiterMax = 500

    while localEff < min_Dtau_frac || localEff > max_Dtau_frac
        if searchiter == 0
            testdistance = 1
        elseif binary_search_mode
            testdistance = (mindist + maxdist) ÷ 2
        else
            testdistance += brute_force_step
        end

        localEff = spatial_decompose(Float64(testdistance),
            dist, tau, E0, nanindex,
            disti, taui, E0i, nanindexi;
            local_eff_mode=true)

        if abs(localEff) < abs(min_Dtau_frac)
            mindist = testdistance
        end
        if abs(localEff) > abs(max_Dtau_frac)
            maxdist = testdistance
        end

        searchiter += 1
        if searchiter == searchiterMax && binary_search_mode
            testdistance = 10
            searchiter = 0
            searchiterMax = Int(round(maximum(filter(!isnan, dist)) / brute_force_step))
            binary_search_mode = false
        elseif searchiter >= searchiterMax
            @warn "Could not find distance solution in find_fraction_distance_localevap"
            return T(NaN)
        end
    end

    return T(testdistance)
end

"""
    find_fraction_distance_streamline(Dtau_bar_val, taubar_init, dx,
        P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
        P0, Fmag0, E0, dist, nanindex, lfrac0;
        dtau_frac=0.75, perc_error=1.0)

Binary search for distance where `dtau_frac` of Dtau_bar is explained
via composite streamlines (upwind effect).
Port of dist_to_Dtaubar_frac_streamline from accessory_fxns.py.
"""
function find_fraction_distance_streamline(Dtau_bar_val::T, taubar_init::T, dx::Real,
        P0i::AbstractVector{T}, Fmag0i::AbstractVector{T},
        E0i::AbstractVector{T}, disti::AbstractVector{T},
        nanindexi::BitVector, lfrac0i::AbstractVector{T},
        P0::AbstractVector{T}, Fmag0::AbstractVector{T},
        E0::AbstractVector{T}, dist::AbstractVector{T},
        nanindex::BitVector, lfrac0::AbstractVector{T};
        dtau_frac::Float64=0.75, perc_error::Float64=1.0) where {T}

    this_Dtau_frac = Dtau_bar_val * dtau_frac
    min_Dtau_frac = minimum([this_Dtau_frac * (1 - perc_error/100),
                             this_Dtau_frac * (1 + perc_error/100)])
    max_Dtau_frac = maximum([this_Dtau_frac * (1 - perc_error/100),
                             this_Dtau_frac * (1 + perc_error/100)])

    mindist = 0
    maxdist = Int(round(maximum(filter(!isnan, dist))))
    localEff = T(-1e3)
    testdistance = 0
    searchiter = 0
    binary_search_mode = true
    brute_force_step = 5
    searchiterMax = 500

    while localEff < min_Dtau_frac || localEff > max_Dtau_frac
        if searchiter == 0
            testdistance = 1
        elseif binary_search_mode
            testdistance = (mindist + maxdist) ÷ 2
        else
            testdistance += brute_force_step
        end

        tau_lc, E0_lc, tau_rc, E0_rc, _, _, _, _, _, _ =
            local_plus_regional_streamline(Float64(testdistance), dx,
                P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
                P0, Fmag0, E0, dist, nanindex, lfrac0)

        taubar_rc = _taubar_simple(tau_rc, E0_rc)
        localEff = taubar_rc - taubar_init

        if abs(localEff) < abs(min_Dtau_frac)
            mindist = testdistance
        end
        if abs(localEff) > abs(max_Dtau_frac)
            maxdist = testdistance
        end

        searchiter += 1
        if searchiter == searchiterMax && binary_search_mode
            testdistance = 10
            searchiter = 0
            searchiterMax = Int(round(maximum(filter(!isnan, dist)) / brute_force_step))
            binary_search_mode = false
        elseif searchiter >= searchiterMax
            @warn "Could not find distance solution in find_fraction_distance_streamline"
            return T(NaN)
        end
    end

    return T(testdistance)
end

# --- E, L, path decomposition ---

"""
    ELpath_decompose(taubar_init,
        E0_E, mu_E, wp_E, tau_E, nanindex_E,
        E0_L, mu_L, wp_L, tau_L, nanindex_L,
        E0_s, mu_s, wp_s, tau_s, nanindex_s)

Decompose tau-bar changes into E (evaporation), L (length scale),
and path (flow field) components.
Port of tau_decompose_ELPath from accessory_fxns.py.

Returns `(Dtau_newE, Dtau_newL, Dtau_newPath)`.
"""
function ELpath_decompose(taubar_init::T,
        E0_E::AbstractVector{T}, mu_E::AbstractVector{T},
        wp_E::AbstractVector{T}, tau_E::AbstractVector{T}, nanindex_E::BitVector,
        E0_L::AbstractVector{T}, mu_L::AbstractVector{T},
        wp_L::AbstractVector{T}, tau_L::AbstractVector{T}, nanindex_L::BitVector,
        E0_s::AbstractVector{T}, mu_s::AbstractVector{T},
        wp_s::AbstractVector{T}, tau_s::AbstractVector{T}, nanindex_s::BitVector) where {T}

    # New E
    mu_tg_E, _, E0_tg_E, ts_E = dist_to_tau_coords(tau_E, nanindex_E, mu_E, wp_E, E0_E)
    taubar_newE = compute_taubar_at_point(mu_tg_E, E0_tg_E, ts_E)
    Dtau_newE = taubar_newE - taubar_init

    # New L
    mu_tg_L, _, E0_tg_L, ts_L = dist_to_tau_coords(tau_L, nanindex_L, mu_L, wp_L, E0_L)
    taubar_newL = compute_taubar_at_point(mu_tg_L, E0_tg_L, ts_L)
    Dtau_newL = taubar_newL - taubar_init

    # New Path
    mu_tg_s, _, E0_tg_s, ts_s = dist_to_tau_coords(tau_s, nanindex_s, mu_s, wp_s, E0_s)
    taubar_news = compute_taubar_at_point(mu_tg_s, E0_tg_s, ts_s)
    Dtau_newPath = taubar_news - taubar_init

    return Dtau_newE, Dtau_newL, Dtau_newPath
end

# --- Main decomposition entry point ---

"""
    decompose_taubar!(ds, config; iso_coords=nothing, runtype="isotopes_noEBM")

Main entry point for tau-bar computation WITH decomposition.
Port of Isotopes.taubar_withDecomposition() from attenuationMod_fxns.py.

Returns `(output_arrays::Dict, df_streamlines::DataFrame)`.
"""
function decompose_taubar!(ds, config::RunConfig;
                           iso_coords::Union{Nothing,DataFrame}=nothing,
                           runtype::String="isotopes_noEBM")
    if occursin("noIsotopes", runtype)
        return Dict{String,Any}(), DataFrame()
    end

    iso = config.isotope
    dc = config.decomp

    # Determine decomposition dimension
    if dc.ELs_initstate_same_yrslice
        decomp_dim = "time"
    else
        decomp_dim = "yr_slice"
    end

    compute_method = iso.compute_method
    streamline_collect = iso.collect_streamline_data
    land_only = iso.compute_land_only
    n_samples = iso.n_samples_taubar_land
    local_threshold = dc.local_threshold_km
    localEvapEff = dc.local_v_regional_local_evap
    UpwindEff = dc.local_v_regional_upwind
    tau_ELs = dc.decomp_E_L_s
    dtau_frac = dc.dtau_fraction

    dx = iso.streamline_dx_km
    Dx = dx / 111.0
    taumax = iso.streamline_max_tau
    dmax = iso.streamline_max_dist_km
    Nmax = ceil(Int, dmax / dx)

    lat = Float64.(ds["lat"][:])
    lon = Float64.(ds["lon"][:])

    # Initialize 3-D arrays
    P, E, UQ, VQ, Fmag, lfrac = tau_array_initialize(ds, lat, lon, config;
                                                       decomp_dim=decomp_dim)

    # Get number of decomposition steps
    n_decomp = size(P, 3)
    baseline_idx = 1  # First index is baseline

    # Build baseline interpolators
    Efiti, Pfiti, uqfiti, vqfiti, Fmagfiti, lfracfiti = build_interpolators(
        lat, lon,
        E[:,:,baseline_idx], P[:,:,baseline_idx],
        UQ[:,:,baseline_idx], VQ[:,:,baseline_idx],
        lfrac[:,:,baseline_idx])

    # Determine computation grid
    if compute_method == "bbox"
        res = iso.bbox_resolution
        comp_lat = collect(Float64, range(-90, 90, step=res[1]))
        comp_lon = collect(Float64, range(0, 360, step=res[2]))
        lat_range = iso.bbox_lat_range
        lon_range = iso.bbox_lon_range
        lat_mask = (comp_lat .>= lat_range[1]) .& (comp_lat .<= lat_range[2])
        lon_mask = (comp_lon .>= lon_range[1]) .& (comp_lon .<= lon_range[2])
        A_slice = findall(lat_mask)
        B_slice = findall(lon_mask)
        grid_lat = comp_lat
        grid_lon = comp_lon
    else
        grid_lat = lat
        grid_lon = lon
        A_slice = Int[]
        B_slice = Int[]
    end

    nlat_g = length(grid_lat)
    nlon_g = length(grid_lon)

    # Output arrays [lat, lon, decomp]
    tau_bar = fill(NaN, nlat_g, nlon_g, n_decomp)
    tau_bar_Ewtd_land = fill(NaN, nlat_g, nlon_g, n_decomp)
    moisture_dist_inland = fill(NaN, nlat_g, nlon_g, n_decomp)
    land_frac_of_streamline = fill(NaN, nlat_g, nlon_g, n_decomp)
    land_frac_of_evapSource = fill(NaN, nlat_g, nlon_g, n_decomp)
    Dtau_bar = fill(NaN, nlat_g, nlon_g, n_decomp)

    Dtau_localEff_localEvap = localEvapEff ? fill(NaN, nlat_g, nlon_g, n_decomp) : nothing
    Dtau_regEff_localEvap = localEvapEff ? fill(NaN, nlat_g, nlon_g, n_decomp) : nothing
    Dtau_localEff_upwind = UpwindEff ? fill(NaN, nlat_g, nlon_g, n_decomp) : nothing
    Dtau_regEff_upwind = UpwindEff ? fill(NaN, nlat_g, nlon_g, n_decomp) : nothing
    Dtau_newE = tau_ELs ? fill(NaN, nlat_g, nlon_g, n_decomp) : nothing
    Dtau_newL = tau_ELs ? fill(NaN, nlat_g, nlon_g, n_decomp) : nothing
    Dtau_newPath = tau_ELs ? fill(NaN, nlat_g, nlon_g, n_decomp) : nothing
    Dtau_50_localEvap = dc.solve_frac_dist_local_evap ? fill(NaN, nlat_g, nlon_g, n_decomp) : nothing
    Dtau_50_strmline = dc.solve_frac_dist_streamline ? fill(NaN, nlat_g, nlon_g, n_decomp) : nothing

    result = StreamlineResult(Nmax)
    resulti = StreamlineResult(Nmax)
    all_streamlines = DataFrame[]

    case_idxs = [i for i in 1:n_decomp if i != baseline_idx]

    if compute_method == "bbox"
        n_cells = length(A_slice) * length(B_slice) * length(case_idxs)
        readout_counter = 1

        for y_idx in A_slice
            for x_idx in B_slice
                lat0 = grid_lat[y_idx]
                lon0 = grid_lon[x_idx]

                land_check = _nearest_landfrac(lfrac[:,:,1], lat, lon, lat0, lon0)
                if land_only && land_check == 0.0
                    moisture_dist_inland[y_idx, x_idx, :] .= 0.0
                    readout_counter += length(case_idxs)
                    continue
                end

                # Baseline streamline
                march_streamline!(resulti, lat0, lon0, Nmax, taumax, Dx, dx,
                                  Efiti, Pfiti, vqfiti, uqfiti, Fmagfiti, lfracfiti)

                tau_wtd_i, land_frac_i, E_terr_i = terrestrial_efrac(
                    land_check, resulti.wp, resulti.coast_step,
                    resulti.tau, resulti.E0, resulti.nanindex, n_samples)

                tau_bar_Ewtd_land[y_idx, x_idx, baseline_idx] = E_terr_i * tau_wtd_i
                max_dist_i = maximum(@view(resulti.dist[1:resulti.nsteps]))
                moisture_dist_inland[y_idx, x_idx, baseline_idx] = resulti.dist_inland
                land_frac_of_streamline[y_idx, x_idx, baseline_idx] = max_dist_i > 0 ? resulti.dist_inland / max_dist_i : 0.0
                land_frac_of_evapSource[y_idx, x_idx, baseline_idx] = land_frac_i

                mu_tg_i, _, E0_tg_i, ts_i = dist_to_tau_coords(
                    resulti.tau, resulti.nanindex, resulti.mu, resulti.wp, resulti.E0)
                tau_bar[y_idx, x_idx, baseline_idx] = compute_taubar_at_point(mu_tg_i, E0_tg_i, ts_i)
                taubar_init = tau_bar[y_idx, x_idx, baseline_idx]

                # Loop through case dimensions
                for thisidx in case_idxs
                    Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = build_interpolators(
                        lat, lon,
                        E[:,:,thisidx], P[:,:,thisidx],
                        UQ[:,:,thisidx], VQ[:,:,thisidx],
                        lfrac[:,:,thisidx])

                    march_streamline!(result, lat0, lon0, Nmax, taumax, Dx, dx,
                                      Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)

                    tau_wtd, land_frac, E_terr = terrestrial_efrac(
                        land_check, result.wp, result.coast_step,
                        result.tau, result.E0, result.nanindex, n_samples)

                    tau_bar_Ewtd_land[y_idx, x_idx, thisidx] = E_terr * tau_wtd
                    max_dist = maximum(@view(result.dist[1:result.nsteps]))
                    moisture_dist_inland[y_idx, x_idx, thisidx] = result.dist_inland
                    land_frac_of_streamline[y_idx, x_idx, thisidx] = max_dist > 0 ? result.dist_inland / max_dist : 0.0
                    land_frac_of_evapSource[y_idx, x_idx, thisidx] = land_frac

                    mu_tg, _, E0_tg, ts = dist_to_tau_coords(
                        result.tau, result.nanindex, result.mu, result.wp, result.E0)
                    tau_bar[y_idx, x_idx, thisidx] = compute_taubar_at_point(mu_tg, E0_tg, ts)
                    Dtau_bar[y_idx, x_idx, thisidx] = tau_bar[y_idx, x_idx, thisidx] - taubar_init

                    # Spatial decomposition: local evap effect
                    if localEvapEff
                        re, le = spatial_decompose(local_threshold,
                            result.dist, result.tau, result.E0, result.nanindex,
                            resulti.dist, resulti.tau, resulti.E0, resulti.nanindex)
                        Dtau_regEff_localEvap[y_idx, x_idx, thisidx] = re
                        Dtau_localEff_localEvap[y_idx, x_idx, thisidx] = le
                    end

                    # Upwind effect
                    if UpwindEff
                        tau_lc, E0_lc, tau_rc, E0_rc, _, _, _, _, _, _ =
                            local_plus_regional_streamline(local_threshold, dx,
                                resulti.P0, resulti.Fmag0, resulti.E0, resulti.dist, resulti.nanindex, resulti.lfrac0,
                                result.P0, result.Fmag0, result.E0, result.dist, result.nanindex, result.lfrac0)
                        taubar_lc = _taubar_simple(tau_lc, E0_lc)
                        taubar_rc = _taubar_simple(tau_rc, E0_rc)
                        Dtau_regEff_upwind[y_idx, x_idx, thisidx] = taubar_lc - taubar_init
                        Dtau_localEff_upwind[y_idx, x_idx, thisidx] = taubar_rc - taubar_init
                    end

                    # E, L, path decomposition
                    if tau_ELs
                        result_E = StreamlineResult(Nmax)
                        result_L = StreamlineResult(Nmax)
                        result_s = StreamlineResult(Nmax)

                        # New E, old L and old flowfield
                        Efit_only, _, _, _, _, _ = build_interpolators(lat, lon,
                            E[:,:,thisidx], P[:,:,baseline_idx],
                            UQ[:,:,baseline_idx], VQ[:,:,baseline_idx], lfrac[:,:,thisidx])
                        march_streamline!(result_E, lat0, lon0, Nmax, taumax, Dx, dx,
                                          Efit_only, Pfiti, vqfiti, uqfiti, Fmagfiti, lfracfit)

                        # New L, old E and old flowfield
                        _, Pfit_only, _, _, Fmagfit_only, _ = build_interpolators(lat, lon,
                            E[:,:,baseline_idx], P[:,:,thisidx],
                            UQ[:,:,baseline_idx], VQ[:,:,baseline_idx], lfrac[:,:,thisidx])
                        march_streamline!(result_L, lat0, lon0, Nmax, taumax, Dx, dx,
                                          Efiti, Pfit_only, vqfiti, uqfiti, Fmagfit_only, lfracfit)

                        # New path, old E and old L
                        _, _, uqfit_only, vqfit_only, _, _ = build_interpolators(lat, lon,
                            E[:,:,baseline_idx], P[:,:,baseline_idx],
                            UQ[:,:,thisidx], VQ[:,:,thisidx], lfrac[:,:,thisidx])
                        march_streamline!(result_s, lat0, lon0, Nmax, taumax, Dx, dx,
                                          Efiti, Pfiti, uqfit_only, vqfit_only, Fmagfiti, lfracfit)

                        dE, dL, ds_val = ELpath_decompose(taubar_init,
                            result_E.E0, result_E.mu, result_E.wp, result_E.tau, result_E.nanindex,
                            result_L.E0, result_L.mu, result_L.wp, result_L.tau, result_L.nanindex,
                            result_s.E0, result_s.mu, result_s.wp, result_s.tau, result_s.nanindex)
                        Dtau_newE[y_idx, x_idx, thisidx] = dE
                        Dtau_newL[y_idx, x_idx, thisidx] = dL
                        Dtau_newPath[y_idx, x_idx, thisidx] = ds_val
                    end

                    # Fraction distance searches
                    if dc.solve_frac_dist_local_evap
                        Dtau_50_localEvap[y_idx, x_idx, thisidx] = find_fraction_distance_localevap(
                            Dtau_bar[y_idx, x_idx, thisidx],
                            result.dist, result.tau, result.E0, result.nanindex,
                            resulti.dist, resulti.tau, resulti.E0, resulti.nanindex;
                            dtau_frac=dtau_frac)
                    end
                    if dc.solve_frac_dist_streamline
                        Dtau_50_strmline[y_idx, x_idx, thisidx] = find_fraction_distance_streamline(
                            Dtau_bar[y_idx, x_idx, thisidx], taubar_init, dx,
                            resulti.P0, resulti.Fmag0, resulti.E0, resulti.dist, resulti.nanindex, resulti.lfrac0,
                            result.P0, result.Fmag0, result.E0, result.dist, result.nanindex, result.lfrac0;
                            dtau_frac=dtau_frac)
                    end

                    pct = round(readout_counter / n_cells * 100, digits=2)
                    println("  decomp tau-bar: $pct% complete")
                    readout_counter += 1
                end
            end
        end

    elseif compute_method == "coord_list" && iso_coords !== nothing
        n_cells = nrow(iso_coords)
        n_iters = n_cells * length(case_idxs)
        readout_counter = 1

        for coord in 1:n_cells
            lat0 = Float64(iso_coords.lat[coord])
            lon0 = Float64(iso_coords.lon[coord])
            if lon0 < 0; lon0 = mod(lon0, 360.0); end

            y_idx = argmin(abs.(lat .- lat0))
            x_idx = argmin(abs.(lon .- lon0))

            land_check = lfrac[y_idx, x_idx, 1]
            if land_only && land_check == 0.0
                moisture_dist_inland[y_idx, x_idx, :] .= 0.0
                readout_counter += length(case_idxs)
                continue
            end

            march_streamline!(resulti, lat0, lon0, Nmax, taumax, Dx, dx,
                              Efiti, Pfiti, vqfiti, uqfiti, Fmagfiti, lfracfiti)

            tau_wtd_i, land_frac_i, E_terr_i = terrestrial_efrac(
                land_check, resulti.wp, resulti.coast_step,
                resulti.tau, resulti.E0, resulti.nanindex, n_samples)

            tau_bar_Ewtd_land[y_idx, x_idx, baseline_idx] = E_terr_i * tau_wtd_i
            max_dist_i = maximum(@view(resulti.dist[1:resulti.nsteps]))
            moisture_dist_inland[y_idx, x_idx, baseline_idx] = resulti.dist_inland
            land_frac_of_streamline[y_idx, x_idx, baseline_idx] = max_dist_i > 0 ? resulti.dist_inland / max_dist_i : 0.0
            land_frac_of_evapSource[y_idx, x_idx, baseline_idx] = land_frac_i

            mu_tg_i, _, E0_tg_i, ts_i = dist_to_tau_coords(
                resulti.tau, resulti.nanindex, resulti.mu, resulti.wp, resulti.E0)
            tau_bar[y_idx, x_idx, baseline_idx] = compute_taubar_at_point(mu_tg_i, E0_tg_i, ts_i)
            taubar_init = tau_bar[y_idx, x_idx, baseline_idx]

            for thisidx in case_idxs
                Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = build_interpolators(
                    lat, lon, E[:,:,thisidx], P[:,:,thisidx],
                    UQ[:,:,thisidx], VQ[:,:,thisidx], lfrac[:,:,thisidx])

                march_streamline!(result, lat0, lon0, Nmax, taumax, Dx, dx,
                                  Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)

                tau_wtd, land_frac, E_terr = terrestrial_efrac(
                    land_check, result.wp, result.coast_step,
                    result.tau, result.E0, result.nanindex, n_samples)

                tau_bar_Ewtd_land[y_idx, x_idx, thisidx] = E_terr * tau_wtd
                max_dist = maximum(@view(result.dist[1:result.nsteps]))
                moisture_dist_inland[y_idx, x_idx, thisidx] = result.dist_inland
                land_frac_of_streamline[y_idx, x_idx, thisidx] = max_dist > 0 ? result.dist_inland / max_dist : 0.0
                land_frac_of_evapSource[y_idx, x_idx, thisidx] = land_frac

                mu_tg, _, E0_tg, ts = dist_to_tau_coords(
                    result.tau, result.nanindex, result.mu, result.wp, result.E0)
                tau_bar[y_idx, x_idx, thisidx] = compute_taubar_at_point(mu_tg, E0_tg, ts)
                Dtau_bar[y_idx, x_idx, thisidx] = tau_bar[y_idx, x_idx, thisidx] - taubar_init

                if localEvapEff
                    re, le = spatial_decompose(local_threshold,
                        result.dist, result.tau, result.E0, result.nanindex,
                        resulti.dist, resulti.tau, resulti.E0, resulti.nanindex)
                    Dtau_regEff_localEvap[y_idx, x_idx, thisidx] = re
                    Dtau_localEff_localEvap[y_idx, x_idx, thisidx] = le
                end

                if UpwindEff
                    tau_lc, E0_lc, tau_rc, E0_rc, _, _, _, _, _, _ =
                        local_plus_regional_streamline(local_threshold, dx,
                            resulti.P0, resulti.Fmag0, resulti.E0, resulti.dist, resulti.nanindex, resulti.lfrac0,
                            result.P0, result.Fmag0, result.E0, result.dist, result.nanindex, result.lfrac0)
                    taubar_lc = _taubar_simple(tau_lc, E0_lc)
                    taubar_rc = _taubar_simple(tau_rc, E0_rc)
                    Dtau_regEff_upwind[y_idx, x_idx, thisidx] = taubar_lc - taubar_init
                    Dtau_localEff_upwind[y_idx, x_idx, thisidx] = taubar_rc - taubar_init
                end

                if tau_ELs
                    result_E = StreamlineResult(Nmax)
                    result_L = StreamlineResult(Nmax)
                    result_s = StreamlineResult(Nmax)

                    Efit_only, _, _, _, _, _ = build_interpolators(lat, lon,
                        E[:,:,thisidx], P[:,:,baseline_idx],
                        UQ[:,:,baseline_idx], VQ[:,:,baseline_idx], lfrac[:,:,thisidx])
                    march_streamline!(result_E, lat0, lon0, Nmax, taumax, Dx, dx,
                                      Efit_only, Pfiti, vqfiti, uqfiti, Fmagfiti, lfracfit)

                    _, Pfit_only, _, _, Fmagfit_only, _ = build_interpolators(lat, lon,
                        E[:,:,baseline_idx], P[:,:,thisidx],
                        UQ[:,:,baseline_idx], VQ[:,:,baseline_idx], lfrac[:,:,thisidx])
                    march_streamline!(result_L, lat0, lon0, Nmax, taumax, Dx, dx,
                                      Efiti, Pfit_only, vqfiti, uqfiti, Fmagfit_only, lfracfit)

                    _, _, uqfit_only, vqfit_only, _, _ = build_interpolators(lat, lon,
                        E[:,:,baseline_idx], P[:,:,baseline_idx],
                        UQ[:,:,thisidx], VQ[:,:,thisidx], lfrac[:,:,thisidx])
                    march_streamline!(result_s, lat0, lon0, Nmax, taumax, Dx, dx,
                                      Efiti, Pfiti, uqfit_only, vqfit_only, Fmagfiti, lfracfit)

                    dE, dL, ds_val = ELpath_decompose(taubar_init,
                        result_E.E0, result_E.mu, result_E.wp, result_E.tau, result_E.nanindex,
                        result_L.E0, result_L.mu, result_L.wp, result_L.tau, result_L.nanindex,
                        result_s.E0, result_s.mu, result_s.wp, result_s.tau, result_s.nanindex)
                    Dtau_newE[y_idx, x_idx, thisidx] = dE
                    Dtau_newL[y_idx, x_idx, thisidx] = dL
                    Dtau_newPath[y_idx, x_idx, thisidx] = ds_val
                end

                pct = round(readout_counter / n_iters * 100, digits=2)
                println("  decomp tau-bar: $pct% complete")
                readout_counter += 1
            end
        end
    end

    df_streamlines = isempty(all_streamlines) ? DataFrame() : vcat(all_streamlines...)

    output = Dict{String,Any}(
        "lat" => grid_lat,
        "lon" => grid_lon,
        "dim3_name" => decomp_dim,
        "tau_bar" => tau_bar,
        "tau_bar_wtdLandEvap" => tau_bar_Ewtd_land,
        "moisture_dist_inland" => moisture_dist_inland,
        "streamline_frac_land" => land_frac_of_streamline,
        "Esource_frac_land" => land_frac_of_evapSource,
        "Dtau_bar" => Dtau_bar,
    )

    if localEvapEff
        output["Dtau_bar_localEff_localEvap"] = Dtau_localEff_localEvap
        output["Dtau_bar_regEff_localEvap"] = Dtau_regEff_localEvap
    end
    if UpwindEff
        output["Dtau_bar_localEff_upwind"] = Dtau_localEff_upwind
        output["Dtau_bar_regEff_upwind"] = Dtau_regEff_upwind
    end
    if tau_ELs
        output["Dtau_bar_Eeffect"] = Dtau_newE
        output["Dtau_bar_Leffect"] = Dtau_newL
        output["Dtau_bar_PathEffect"] = Dtau_newPath
    end
    if dc.solve_frac_dist_local_evap
        output["Dtau_bar_fracDist_localEvap"] = Dtau_50_localEvap
    end
    if dc.solve_frac_dist_streamline
        output["Dtau_bar_fracDist_upwind"] = Dtau_50_strmline
    end

    return output, df_streamlines
end
