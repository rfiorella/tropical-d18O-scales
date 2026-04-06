using KernelAbstractions
using StaticArrays

"""
GPU-compatible bilinear interpolation using precomputed grid parameters from DeltapFields.
"""
@inline function bilinear_interp_gpu(field, fields::DeltapFields{T}, lat, lon) where T
    fi = (lat - fields.lat0) * fields.inv_dlat
    fj = (lon - fields.lon0) * fields.inv_dlon

    i = clamp(unsafe_trunc(Int32, fi) + Int32(1), Int32(1), fields.n_lat - Int32(1))
    j = clamp(unsafe_trunc(Int32, fj) + Int32(1), Int32(1), fields.n_lon - Int32(1))

    t_lat = clamp(fi - T(i - Int32(1)), T(0), T(1))
    t_lon = clamp(fj - T(j - Int32(1)), T(0), T(1))

    @inbounds begin
        return (one(T) - t_lat) * (one(T) - t_lon) * field[i, j] +
               (one(T) - t_lat) *           t_lon  * field[i, j+1] +
                         t_lat  * (one(T) - t_lon) * field[i+1, j] +
                         t_lat  *           t_lon  * field[i+1, j+1]
    end
end

"""
GPU-compatible 1D interpolation for alpha_eq lookup on regularly-spaced table.
"""
@inline function alpha_lookup_gpu(fields::DeltapFields{T}, temp_c) where T
    fi = (temp_c - fields.alpha_Tmin) * fields.alpha_inv_dT
    n = Int32(length(fields.alpha_T))
    lo = clamp(unsafe_trunc(Int32, fi) + Int32(1), Int32(1), n - Int32(1))
    t = fi - T(lo - Int32(1))
    @inbounds return fields.alpha_val[lo] + t * (fields.alpha_val[lo+1] - fields.alpha_val[lo])
end

"""
Walk a streamline from (lat0, lon0), returning (lat_new, lon_new, vq, uq) at each step.
Inline helper to avoid duplicating the advection logic between passes.
"""
@inline function advect_step(fields::DeltapFields{T}, lat_cur, lon_cur,
                              vq_cur, uq_cur, fm_cur, Dx) where T
    inv_fm = T(1) / fm_cur
    dtheta = -vq_cur * inv_fm * Dx
    dphi   = -uq_cur * inv_fm * Dx / cos((lat_cur + dtheta / T(2)) * T(π) / T(180))

    lat_new = lat_cur + dtheta
    lon_new = lon_cur + dphi

    if lat_new < T(-90)
        lat_new = T(180) + lat_new; lon_new -= T(180)
    elseif lat_new > T(90)
        lat_new = T(180) - lat_new; lon_new -= T(180)
    end
    if lon_new < T(0)
        lon_new += T(360)
    elseif lon_new > T(360)
        lon_new -= T(360)
    end

    return lat_new, lon_new
end

"""
Two-trace GPU kernel: re-walks the streamline twice to eliminate all per-step array storage.
Pass 1: trace streamline, compute npts and total_mu_sum (for decay weights).
Pass 2: re-trace, compute weighted averages in a single forward pass.

This uses O(1) registers per thread instead of O(Nmax), eliminating register spill on Metal.
"""
@kernel function deltap_kernel!(deltap_bar, tau_bar,
                                 @Const(LAT2), @Const(LON2), @Const(P2), @Const(Plim2),
                                 fields::DeltapFields{T},
                                 dx::T, Dx::T, Nmax_val::Int32, tau_stop::T) where T

    idx = @index(Global, Linear)
    A = Int32(size(LAT2, 1))
    j_out = (idx - Int32(1)) ÷ A + Int32(1)
    i_out = (idx - Int32(1)) % A + Int32(1)

    @inbounds if P2[i_out, j_out] >= Plim2[i_out, j_out]

    lat_start = T(@inbounds LAT2[i_out, j_out])
    lon_start = T(@inbounds LON2[i_out, j_out])

    dx1000 = dx * T(1000)

    # =====================================================
    # PASS 1: Trace streamline, compute npts + total_mu_sum
    # =====================================================
    lat_cur = lat_start
    lon_cur = lon_start

    vq_cur = bilinear_interp_gpu(fields.VQ, fields, lat_cur, lon_cur)
    uq_cur = bilinear_interp_gpu(fields.UQ, fields, lat_cur, lon_cur)
    fm_cur = sqrt(uq_cur^2 + vq_cur^2)
    p_cur  = bilinear_interp_gpu(fields.P, fields, lat_cur, lon_cur)
    de_cur = bilinear_interp_gpu(fields.delta_e, fields, lat_cur, lon_cur)

    # tau and mu for step 1
    mu_cur = p_cur / fm_cur
    tau_running = mu_cur * dx1000  # tau[1] after recomputation
    tau_prev = T(0)  # for the while-loop tau (different from post-processing tau)

    # total_mu_sum accumulation (step 1 is always valid)
    total_mu_sum = mu_cur * dx1000

    npts = Int32(1)
    jj = Int32(1)

    while tau_prev < tau_stop && jj < Nmax_val - Int32(1)
        lat_new, lon_new = advect_step(fields, lat_cur, lon_cur, vq_cur, uq_cur, fm_cur, Dx)

        vq_new = bilinear_interp_gpu(fields.VQ, fields, lat_new, lon_new)
        uq_new = bilinear_interp_gpu(fields.UQ, fields, lat_new, lon_new)
        fm_new = sqrt(uq_new^2 + vq_new^2)
        p_new  = bilinear_interp_gpu(fields.P, fields, lat_new, lon_new)
        de_new = bilinear_interp_gpu(fields.delta_e, fields, lat_new, lon_new)

        # Advance the while-loop tau (uses previous step's P/F)
        mu_val = p_cur / fm_cur
        tau_prev = tau_prev + mu_val * dx1000

        jj += Int32(1)
        npts = jj

        # Post-processing tau for this step
        mu_new = p_new / fm_new
        tau_running += mu_new * dx1000  # cumulative tau up to step jj

        # Check validity for mu_sum
        tau_k = tau_running  # this is the post-processing tau at step jj
        is_valid = tau_k > T(0) && de_new < T(120) && de_new > T(-120)
        if is_valid
            total_mu_sum += mu_new * dx1000
        end

        vq_cur = vq_new; uq_cur = uq_new; fm_cur = fm_new
        p_cur = p_new; de_cur = de_new
        lat_cur = lat_new; lon_cur = lon_new
    end

    # =====================================================
    # PASS 2: Re-trace, compute weighted delta_p and tau
    # =====================================================
    lat_cur = lat_start
    lon_cur = lon_start

    vq_cur = bilinear_interp_gpu(fields.VQ, fields, lat_cur, lon_cur)
    uq_cur = bilinear_interp_gpu(fields.UQ, fields, lat_cur, lon_cur)
    fm_cur = sqrt(uq_cur^2 + vq_cur^2)
    p_cur  = bilinear_interp_gpu(fields.P, fields, lat_cur, lon_cur)
    tc_cur = bilinear_interp_gpu(fields.Tcond, fields, lat_cur, lon_cur)
    e_cur  = bilinear_interp_gpu(fields.E, fields, lat_cur, lon_cur)
    de_cur = bilinear_interp_gpu(fields.delta_e, fields, lat_cur, lon_cur)

    # Epsilon at origin
    epsilonc0_1 = (alpha_lookup_gpu(fields, tc_cur - T(273.15)) - T(1)) * T(1000)

    # Step 1 values
    mu_k = p_cur / fm_cur
    tau_k = mu_k * dx1000  # tau[1]
    mu_dx_k = mu_k * dx1000

    # Step 1 is always valid → start accumulation
    cum_mu = T(0)
    decay_1 = exp(-(total_mu_sum - cum_mu))
    cum_mu += mu_dx_k
    w_1 = mu_k * decay_1

    alpha_1 = alpha_lookup_gpu(fields, tc_cur - T(273.15))
    cum_weight_alpha = w_1 * alpha_1
    cum_weight = w_1

    alpha_bar_1 = cum_weight_alpha / cum_weight
    epsilon_bar_1 = (alpha_bar_1 - T(1)) * T(1000)
    delta_p_1 = -tau_k * epsilon_bar_1 + epsilonc0_1 + de_cur

    exp_neg_tau_1 = exp(-tau_k)
    ew_1 = e_cur * exp_neg_tau_1

    num_dp = delta_p_1 * ew_1
    num_t = tau_k * ew_1
    den = ew_1

    # Walk remaining steps
    tau_running2 = tau_k
    for step in Int32(2):npts
        lat_new, lon_new = advect_step(fields, lat_cur, lon_cur, vq_cur, uq_cur, fm_cur, Dx)

        vq_new = bilinear_interp_gpu(fields.VQ, fields, lat_new, lon_new)
        uq_new = bilinear_interp_gpu(fields.UQ, fields, lat_new, lon_new)
        fm_new = sqrt(uq_new^2 + vq_new^2)
        p_new  = bilinear_interp_gpu(fields.P, fields, lat_new, lon_new)
        tc_new = bilinear_interp_gpu(fields.Tcond, fields, lat_new, lon_new)
        e_new  = bilinear_interp_gpu(fields.E, fields, lat_new, lon_new)
        de_new = bilinear_interp_gpu(fields.delta_e, fields, lat_new, lon_new)

        mu_k = p_new / fm_new
        mu_dx_k = mu_k * dx1000
        tau_running2 += mu_dx_k
        tau_k = tau_running2

        is_valid = tau_k > T(0) && de_new < T(120) && de_new > T(-120)
        if is_valid
            decay_k = exp(-(total_mu_sum - cum_mu))
            cum_mu += mu_dx_k
            w_k = mu_k * decay_k

            alpha_k = alpha_lookup_gpu(fields, tc_new - T(273.15))
            cum_weight_alpha += w_k * alpha_k
            cum_weight += w_k

            alpha_bar_k = cum_weight > T(0) ? cum_weight_alpha / cum_weight : T(1.01)
            epsilon_bar_k = (alpha_bar_k - T(1)) * T(1000)
            delta_p_k = -tau_k * epsilon_bar_k + epsilonc0_1 + de_new

            exp_neg_tau_k = exp(-tau_k)
            ew_k = e_new * exp_neg_tau_k

            num_dp += delta_p_k * ew_k
            num_t += tau_k * ew_k
            den += ew_k
        end

        vq_cur = vq_new; uq_cur = uq_new; fm_cur = fm_new
        p_cur = p_new
        lat_cur = lat_new; lon_cur = lon_new
    end

    if den > T(0)
        @inbounds deltap_bar[i_out, j_out] = num_dp / den
        @inbounds tau_bar[i_out, j_out] = num_t / den
    end

    end # if P2 >= Plim2
end
