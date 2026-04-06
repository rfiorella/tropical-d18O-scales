"""
    select_backend()

Auto-select GPU backend: Metal → CPU fallback.
Returns (backend, ArrayType).
"""
function select_backend()
    try
        if Metal.functional()
            return Metal.MetalBackend(), Metal.MtlArray
        end
    catch
    end
    return CPU(), Array
end

"""
    run_gpu(E, P, UQ, VQ, Tcond, LAT, LON, LAT2, LON2,
            delta_e, alpha_eq, Plim, dmax;
            T=Float32, workgroup_size=256, backend=nothing)

GPU-accelerated version of get_deltap_cpu.

Prepares data, uploads to GPU, launches deltap_kernel!, downloads results.
Falls back to KernelAbstractions CPU backend if no GPU is available.

# Keyword arguments
- `T`: Float type (default Float32 for GPU bandwidth)
- `workgroup_size`: GPU workgroup size (default 256)
- `backend`: Override backend (default: auto-select)
"""
function run_gpu(E, P, UQ, VQ, Tcond, LAT, LON, LAT2, LON2,
                 delta_e, alpha_eq, Plim, dmax;
                 T=Float32, workgroup_size=256, backend=nothing, fixed_nsteps::Int = 0)

    if backend === nothing
        be, ArrayType = select_backend()
    else
        be = backend
        ArrayType = be isa KernelAbstractions.CPU ? Array : Metal.MtlArray
    end

    #---------------------------------------
    # Setup (same as CPU version)
    #---------------------------------------
    dx  = T(15.0)
    Dx  = T(dx / 111.0)
    Nmax = fixed_nsteps > 0 ? fixed_nsteps : ceil(Int, dmax / 15.0)
    tau_stop = T(fixed_nsteps > 0 ? Inf : 10.0)

    # Clamp and inpaint
    delta_e = copy(delta_e)
    delta_e .= clamp.(delta_e, -200.0, 200.0)
    Tcond = copy(Tcond)
    if any(isnan, delta_e)
        delta_e = Float64.(inpaint_nans!(Float64.(delta_e)))
    end
    if any(isnan, Tcond)
        Tcond = Float64.(inpaint_nans!(Float64.(Tcond)))
    end

    # Cyclic extension along dim 1 (longitude)
    LON_ext = vcat(LON[end:end, :] .- 360, LON, LON[1:1, :] .+ 360)
    LAT_ext = vcat(LAT[end:end, :], LAT, LAT[1:1, :])
    E_ext   = vcat(E[end:end, :], E, E[1:1, :])
    P_ext   = vcat(P[end:end, :], P, P[1:1, :])
    UQ_ext  = vcat(UQ[end:end, :], UQ, UQ[1:1, :])
    VQ_ext  = vcat(VQ[end:end, :], VQ, VQ[1:1, :])
    Tcond_ext = vcat(Tcond[end:end, :], Tcond, Tcond[1:1, :])
    Plim_ext  = vcat(Plim[end:end, :], Plim, Plim[1:1, :])
    delta_e_ext = vcat(delta_e[end:end, :], delta_e, delta_e[1:1, :])

    # Transpose to (nlat, nlon) for interpolation
    P_interp       = T.(permutedims(P_ext))
    UQ_interp      = T.(permutedims(UQ_ext))
    VQ_interp      = T.(permutedims(VQ_ext))
    E_interp       = T.(permutedims(E_ext))
    Tcond_interp   = T.(permutedims(Tcond_ext))
    delta_e_interp = T.(permutedims(delta_e_ext))
    Plim_interp    = T.(permutedims(Plim_ext))

    lat_vec = T.(LAT_ext[1, :])
    lon_vec = T.(LON_ext[:, 1])

    dlat = lat_vec[2] - lat_vec[1]
    dlon = lon_vec[2] - lon_vec[1]

    # Alpha_eq precomputation
    sort_idx = sortperm(alpha_eq[:, 1])
    alpha_T_sorted   = T.(alpha_eq[sort_idx, 1])
    alpha_val_sorted = T.(alpha_eq[sort_idx, 2])
    alpha_dT = alpha_T_sorted[2] - alpha_T_sorted[1]

    # Build DeltapFields struct
    fields = DeltapFields(
        ArrayType(P_interp), ArrayType(UQ_interp), ArrayType(VQ_interp),
        ArrayType(E_interp), ArrayType(Tcond_interp),
        ArrayType(delta_e_interp), ArrayType(Plim_interp),
        T(lat_vec[1]), T(lon_vec[1]),
        T(1.0 / dlat), T(1.0 / dlon),
        Int32(length(lat_vec)), Int32(length(lon_vec)),
        ArrayType(alpha_T_sorted), ArrayType(alpha_val_sorted),
        T(alpha_T_sorted[1]), T(1.0 / alpha_dT)
    )

    # Precompute P2 and Plim2 on CPU (small cost), then upload
    grid = RegularGrid(Float64.(lat_vec), Float64.(lon_vec))
    A, B = size(LAT2)
    P2_cpu    = zeros(T, A, B)
    Plim2_cpu = zeros(T, A, B)
    for jj in 1:B, ii in 1:A
        P2_cpu[ii, jj]    = T(bilinear_interp(Float64.(P_interp), grid, Float64(LAT2[ii, jj]), Float64(LON2[ii, jj])))
        Plim2_cpu[ii, jj] = T(bilinear_interp(Float64.(Plim_interp), grid, Float64(LAT2[ii, jj]), Float64(LON2[ii, jj])))
    end

    # Upload to device
    d_LAT2 = ArrayType(T.(LAT2))
    d_LON2 = ArrayType(T.(LON2))
    d_P2   = ArrayType(P2_cpu)
    d_Plim2 = ArrayType(Plim2_cpu)
    d_deltap = ArrayType(fill(T(NaN), A, B))
    d_tau    = ArrayType(fill(T(NaN), A, B))

    # Launch kernel
    kernel! = deltap_kernel!(be, workgroup_size)
    kernel!(d_deltap, d_tau, d_LAT2, d_LON2, d_P2, d_Plim2,
            fields, dx, Dx, Int32(Nmax), tau_stop;
            ndrange=A * B)
    KernelAbstractions.synchronize(be)

    # Download and convert to Float64 on CPU (Metal doesn't support Float64)
    return Float64.(Array(d_deltap)), Float64.(Array(d_tau))
end
