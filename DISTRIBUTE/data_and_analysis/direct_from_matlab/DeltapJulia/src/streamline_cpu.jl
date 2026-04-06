"""
    get_deltap_cpu(E, P, UQ, VQ, Tcond, LAT, LON, LAT2, LON2,
                   delta_e, alpha_eq, Plim, dmax)

Lagrangian streamline-based Rayleigh distillation model for precipitation δ¹⁸O.
Line-for-line translation of get_deltap_fast_optimized.m.

# Input convention (matching Matlab)
All 2D fields have size (nlon, nlat):
- Dimension 1 (rows) = longitude
- Dimension 2 (cols) = latitude

`alpha_eq` is an Nx2 matrix: col 1 = temperature (°C), col 2 = equilibrium α.

# Returns
- `deltap_bar`: precipitation δ¹⁸O (same size as LAT2)
- `tau_bar`: residence-time metric (same size as LAT2)
"""
function get_deltap_cpu(E, P, UQ, VQ, Tcond, LAT, LON, LAT2, LON2,
                        delta_e, alpha_eq, Plim, dmax;
                        fixed_nsteps::Int = 0)

    #---------------------------------------
    # Setup constants and grid extensions
    #---------------------------------------
    dx  = 15.0          # km
    Dx  = dx / 111.0    # deg latitude increment
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
    # Matlab: cat(1, X(end,:), X, X(1,:)) — wraps longitude
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
    # Build interpolation grids
    #---------------------------------------
    # Matlab does griddedInterpolant(LAT', LON', F') which transposes
    # from (nlon, nlat) -> (nlat, nlon), making lat the row dimension.
    # We transpose to (nlat, nlon) for Julia bilinear_interp:
    #   lat_vec = sorted latitudes (row dimension)
    #   lon_vec = sorted longitudes (col dimension)
    P_interp       = permutedims(P_ext)
    UQ_interp      = permutedims(UQ_ext)
    VQ_interp      = permutedims(VQ_ext)
    E_interp       = permutedims(E_ext)
    Tcond_interp   = permutedims(Tcond_ext)
    delta_e_interp = permutedims(delta_e_ext)
    Plim_interp    = permutedims(Plim_ext)

    # 1D coordinate vectors from the transposed grid
    lat_vec = LAT_ext[1, :]   # latitude varies along original cols = new rows after transpose
    lon_vec = LON_ext[:, 1]   # longitude varies along original rows = new cols after transpose

    # After transpose: field[i,j] where i indexes lat_vec, j indexes lon_vec
    # So we need lat_vec as the "row" coords and lon_vec as the "col" coords
    # But wait — permutedims swaps dims, so:
    #   P_interp[i, j] = P_ext[j, i]
    #   i ranges over 1:nlat, j ranges over 1:nlon_ext
    # lat_vec should be LAT_ext[1, :] = latitudes (size nlat)
    # lon_vec should be LON_ext[:, 1] = longitudes (size nlon_ext)

    # Helper: interpolate a field at (lat, lon)
    function do_interp(field, lat, lon)
        bilinear_interp(field, lat_vec, lon_vec, lat, lon)
    end

    #---------------------------------------
    # Initialize outputs
    #---------------------------------------
    A, B = size(LAT2)
    deltap_bar = fill(NaN, A, B)
    tau_bar    = fill(NaN, A, B)

    #---------------------------------------
    # Precompute P2 and Plim2 for validity check
    #---------------------------------------
    P2    = similar(LAT2, Float64)
    Plim2 = similar(LAT2, Float64)
    for j in 1:B, i in 1:A
        P2[i, j]    = do_interp(P_interp, LAT2[i, j], LON2[i, j])
        Plim2[i, j] = do_interp(Plim_interp, LAT2[i, j], LON2[i, j])
    end

    # Sort alpha_eq by temperature (ascending) for interp1
    sort_idx = sortperm(alpha_eq[:, 1])
    alpha_T   = alpha_eq[sort_idx, 1]
    alpha_val = alpha_eq[sort_idx, 2]

    #---------------------------------------
    # Main grid loops
    #---------------------------------------
    for i in 1:A
        for j in 1:B

            # Skip invalid points
            if P2[i, j] < Plim2[i, j]
                continue
            end

            #---------------------------------------
            # Initialize streamline integration
            #---------------------------------------
            lat0 = Float64(LAT2[i, j])
            lon0 = Float64(LON2[i, j])

            tau      = zeros(Nmax)
            E0       = zeros(Nmax)
            Tcond0   = zeros(Nmax)
            Fmag0    = zeros(Nmax)
            P0       = zeros(Nmax)
            delta_e0 = zeros(Nmax)

            vq0_val = do_interp(VQ_interp, lat0, lon0)
            uq0_val = do_interp(UQ_interp, lat0, lon0)
            E0[1]       = do_interp(E_interp, lat0, lon0)
            Tcond0[1]   = do_interp(Tcond_interp, lat0, lon0)
            Fmag0[1]    = hypot(uq0_val, vq0_val)
            P0[1]       = do_interp(P_interp, lat0, lon0)
            delta_e0[1] = do_interp(delta_e_interp, lat0, lon0)

            jj = 1
            while tau[jj] < tau_stop && jj < Nmax - 1
                # Streamline step
                dtheta = -vq0_val / Fmag0[jj] * Dx
                dphi   = -uq0_val / Fmag0[jj] * Dx / cosd(lat0 + dtheta / 2)

                lat1 = lat0 + dtheta
                lon1 = lon0 + dphi

                # Latitude/longitude wrapping
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

                # Interpolate fields at new point
                vq1 = do_interp(VQ_interp, lat1, lon1)
                uq1 = do_interp(UQ_interp, lat1, lon1)
                E0[jj+1]       = do_interp(E_interp, lat1, lon1)
                Tcond0[jj+1]   = do_interp(Tcond_interp, lat1, lon1)
                Fmag0[jj+1]    = hypot(uq1, vq1)
                P0[jj+1]       = do_interp(P_interp, lat1, lon1)
                delta_e0[jj+1] = do_interp(delta_e_interp, lat1, lon1)

                # Advance
                vq0_val = vq1; uq0_val = uq1
                lat0 = lat1; lon0 = lon1

                jj += 1
                mu_val = P0[jj-1] / Fmag0[jj-1]
                tau[jj] = tau[jj-1] + mu_val * dx * 1000
            end

            #---------------------------------------
            # Compute along-streamline properties
            #---------------------------------------
            mu = P0 ./ Fmag0
            tau = cumsum(mu) .* dx * 1000

            valid_idx = (tau .> 0) .& (delta_e0 .< 120) .& (delta_e0 .> -120)
            valid_idx[1] = true

            # Alpha_eq lookup via 1D interpolation with extrapolation
            alphac0 = interp1_extrap(alpha_T, alpha_val, Tcond0 .- 273.15)
            epsilonc0 = (alphac0 .- 1) .* 1000

            # Weight function
            weight_core = mu[valid_idx]
            decay = exp.(-reverse(cumsum(reverse(weight_core) .* dx * 1000)))
            weight = vcat(weight_core .* decay, zeros(sum(.!valid_idx)))

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
