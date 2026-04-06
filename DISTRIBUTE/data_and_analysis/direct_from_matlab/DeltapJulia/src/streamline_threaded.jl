"""
    get_deltap_cpu_threaded(E, P, UQ, VQ, Tcond, LAT, LON, LAT2, LON2,
                            delta_e, alpha_eq, Plim, dmax)

Threaded, optimized version of `get_deltap_cpu`.
Uses O(1) regular-grid interpolation, precomputed alpha_eq lookup,
`Threads.@threads` on the outer loop, and pre-allocated per-thread workspaces.

Same inputs/outputs as `get_deltap_cpu`.
"""
function get_deltap_cpu_threaded(E, P, UQ, VQ, Tcond, LAT, LON, LAT2, LON2,
                                  delta_e, alpha_eq, Plim, dmax;
                                  fixed_nsteps::Int = 0)

    #---------------------------------------
    # Setup constants and grid extensions
    #---------------------------------------
    dx  = 15.0
    Dx  = dx / 111.0
    Nmax = fixed_nsteps > 0 ? fixed_nsteps : ceil(Int, dmax / dx)
    tau_stop = fixed_nsteps > 0 ? Inf : 10.0

    # Clamp delta_e
    delta_e = copy(delta_e)
    delta_e .= clamp.(delta_e, -200.0, 200.0)

    # Remove NaNs from delta_e and Tcond
    Tcond = copy(Tcond)
    if any(isnan, delta_e)
        delta_e = Float64.(inpaint_nans!(Float64.(delta_e)))
    end
    if any(isnan, Tcond)
        Tcond = Float64.(inpaint_nans!(Float64.(Tcond)))
    end

    # Expand grid cyclically along dimension 1 (longitude rows)
    LON_ext = vcat(LON[end:end, :] .- 360, LON, LON[1:1, :] .+ 360)
    LAT_ext = vcat(LAT[end:end, :], LAT, LAT[1:1, :])
    E_ext   = vcat(E[end:end, :], E, E[1:1, :])
    P_ext   = vcat(P[end:end, :], P, P[1:1, :])
    UQ_ext  = vcat(UQ[end:end, :], UQ, UQ[1:1, :])
    VQ_ext  = vcat(VQ[end:end, :], VQ, VQ[1:1, :])
    Tcond_ext = vcat(Tcond[end:end, :], Tcond, Tcond[1:1, :])
    Plim_ext  = vcat(Plim[end:end, :], Plim, Plim[1:1, :])
    delta_e_ext = vcat(delta_e[end:end, :], delta_e, delta_e[1:1, :])

    #---------------------------------------
    # Build interpolation grids (transposed)
    #---------------------------------------
    P_interp       = permutedims(P_ext)
    UQ_interp      = permutedims(UQ_ext)
    VQ_interp      = permutedims(VQ_ext)
    E_interp       = permutedims(E_ext)
    Tcond_interp   = permutedims(Tcond_ext)
    delta_e_interp = permutedims(delta_e_ext)
    Plim_interp    = permutedims(Plim_ext)

    lat_vec = LAT_ext[1, :]
    lon_vec = LON_ext[:, 1]

    # Precompute RegularGrid for O(1) index lookup
    grid = RegularGrid(lat_vec, lon_vec)

    #---------------------------------------
    # Precompute alpha_eq lookup
    #---------------------------------------
    sort_idx = sortperm(alpha_eq[:, 1])
    alpha_T   = alpha_eq[sort_idx, 1]
    alpha_val = alpha_eq[sort_idx, 2]
    alpha_lut = RegularLookup(alpha_T, alpha_val)

    #---------------------------------------
    # Precompute P2 and Plim2
    #---------------------------------------
    A, B = size(LAT2)
    P2    = similar(LAT2, Float64)
    Plim2 = similar(LAT2, Float64)
    for j in 1:B, i in 1:A
        P2[i, j]    = bilinear_interp(P_interp, grid, LAT2[i, j], LON2[i, j])
        Plim2[i, j] = bilinear_interp(Plim_interp, grid, LAT2[i, j], LON2[i, j])
    end

    #---------------------------------------
    # Initialize outputs
    #---------------------------------------
    deltap_bar = fill(NaN, A, B)
    tau_bar    = fill(NaN, A, B)

    #---------------------------------------
    # Per-thread workspace allocation
    #---------------------------------------
    nthreads = max(Threads.nthreads(), Threads.maxthreadid())
    ws_tau      = [zeros(Nmax) for _ in 1:nthreads]
    ws_E0       = [zeros(Nmax) for _ in 1:nthreads]
    ws_Tcond0   = [zeros(Nmax) for _ in 1:nthreads]
    ws_Fmag0    = [zeros(Nmax) for _ in 1:nthreads]
    ws_P0       = [zeros(Nmax) for _ in 1:nthreads]
    ws_delta_e0 = [zeros(Nmax) for _ in 1:nthreads]

    #---------------------------------------
    # Threaded main loop (column-major: outer=j)
    #---------------------------------------
    Threads.@threads for j in 1:B
        tid = Threads.threadid()
        tau      = ws_tau[tid]
        E0       = ws_E0[tid]
        Tcond0   = ws_Tcond0[tid]
        Fmag0    = ws_Fmag0[tid]
        P0       = ws_P0[tid]
        delta_e0 = ws_delta_e0[tid]

        for i in 1:A
            if P2[i, j] < Plim2[i, j]
                continue
            end

            #---------------------------------------
            # Initialize streamline integration
            #---------------------------------------
            lat0 = Float64(LAT2[i, j])
            lon0 = Float64(LON2[i, j])

            fill!(tau, 0.0)
            fill!(E0, 0.0)
            fill!(Tcond0, 0.0)
            fill!(Fmag0, 0.0)
            fill!(P0, 0.0)
            fill!(delta_e0, 0.0)

            vq0_val = bilinear_interp(VQ_interp, grid, lat0, lon0)
            uq0_val = bilinear_interp(UQ_interp, grid, lat0, lon0)
            E0[1]       = bilinear_interp(E_interp, grid, lat0, lon0)
            Tcond0[1]   = bilinear_interp(Tcond_interp, grid, lat0, lon0)
            Fmag0[1]    = hypot(uq0_val, vq0_val)
            P0[1]       = bilinear_interp(P_interp, grid, lat0, lon0)
            delta_e0[1] = bilinear_interp(delta_e_interp, grid, lat0, lon0)

            jj = 1
            while tau[jj] < tau_stop && jj < Nmax - 1
                dtheta = -vq0_val / Fmag0[jj] * Dx
                dphi   = -uq0_val / Fmag0[jj] * Dx / cosd(lat0 + dtheta / 2)

                lat1 = lat0 + dtheta
                lon1 = lon0 + dphi

                if lat1 < -90
                    lat1 = 180 + lat1; lon1 = lon1 - 180
                elseif lat1 > 90
                    lat1 = 180 - lat1; lon1 = lon1 - 180
                end
                if lon1 < 0
                    lon1 = lon1 + 360
                elseif lon1 > 360
                    lon1 = lon1 - 360
                end

                vq1 = bilinear_interp(VQ_interp, grid, lat1, lon1)
                uq1 = bilinear_interp(UQ_interp, grid, lat1, lon1)
                E0[jj+1]       = bilinear_interp(E_interp, grid, lat1, lon1)
                Tcond0[jj+1]   = bilinear_interp(Tcond_interp, grid, lat1, lon1)
                Fmag0[jj+1]    = hypot(uq1, vq1)
                P0[jj+1]       = bilinear_interp(P_interp, grid, lat1, lon1)
                delta_e0[jj+1] = bilinear_interp(delta_e_interp, grid, lat1, lon1)

                vq0_val = vq1; uq0_val = uq1
                lat0 = lat1; lon0 = lon1

                jj += 1
                mu_val = P0[jj-1] / Fmag0[jj-1]
                tau[jj] = tau[jj-1] + mu_val * dx * 1000
            end

            #---------------------------------------
            # Post-processing (vectorized, matches naive version)
            #---------------------------------------
            mu = P0 ./ Fmag0
            @. tau = 0.0
            tau[1] = mu[1] * dx * 1000
            for k in 2:Nmax
                tau[k] = tau[k-1] + mu[k] * dx * 1000
            end

            valid_idx = (tau .> 0) .& (delta_e0 .< 120) .& (delta_e0 .> -120)
            valid_idx[1] = true

            # Alpha_eq lookup via precomputed regular lookup
            alphac0 = [interp1_extrap(alpha_lut, Tcond0[k] - 273.15) for k in 1:Nmax]
            epsilonc0 = (alphac0 .- 1) .* 1000

            # Weight function
            weight_core = mu[valid_idx]
            decay = exp.(-reverse(cumsum(reverse(weight_core) .* dx * 1000)))
            weight = vcat(weight_core .* decay, zeros(count(.!valid_idx)))

            # Mean alpha and epsilon
            alpha_bar0 = cumsum(weight .* alphac0) ./ cumsum(weight)
            alpha_bar0[isnan.(alpha_bar0)] .= 1.01
            epsilon_bar0 = (alpha_bar0 .- 1) .* 1000

            # delta_p and exponential weighting
            delta_p0 = -tau .* epsilon_bar0 .+ epsilonc0[1] .+ delta_e0
            exp_neg_tau = exp.(-tau[valid_idx])
            eweight = E0[valid_idx] .* exp_neg_tau

            num_deltap = sum(delta_p0[valid_idx] .* eweight)
            num_tau    = sum(tau[valid_idx] .* eweight)
            den        = sum(eweight)

            if den > 0
                deltap_bar[i, j] = num_deltap / den
                tau_bar[i, j]    = num_tau / den
            end
        end
    end

    return deltap_bar, tau_bar
end
