# Task Spec: Matlab → Julia GPU Port of `get_deltap_fast_optimized`

## Overview

Port a Lagrangian streamline-based Rayleigh distillation model for precipitation δ¹⁸O
from Matlab to Julia, with GPU acceleration targeting both NVIDIA (CUDA) and Apple
Silicon (Metal) from a single codebase.

The Matlab function traces moisture backward along integrated vapor transport (IVT)
streamlines for each output grid point, accumulating isotopic fractionation and mixing
with evaporative source vapor to compute precipitation δ¹⁸O (`deltap_bar`) and a
residence-time metric (`tau_bar`).

**Source file**: `get_deltap_fast_optimized.m` (included below in Appendix A)

---

## Computational Profile of the Source Code

### Structure

```
for i = 1:A          ← outer grid loop (embarrassingly parallel)
  for j = 1:B        ← outer grid loop (embarrassingly parallel)
    if P2(i,j) < Plim2(i,j): skip

    while tau < 10 && jj < Nmax:   ← streamline tracing (serial per point)
      - bilinear interpolation of 6 fields (UQ, VQ, E, P, Tcond, delta_e)
      - advection step (lat/lon update with wrapping)
      - tau accumulation

    post-processing:               ← per-point reduction
      - cumsum for tau
      - 1D interp1 for alpha_eq lookup
      - exponential-decay weighting
      - weighted mean for deltap_bar, tau_bar
```

### Key characteristics

- **Parallelism**: The `(i,j)` outer loop is fully independent — no data dependencies
  between grid points. This is the GPU parallelization target.
- **Inner loop**: Variable-length `while` loop with ceiling at `Nmax = ceil(dmax/15)`.
  Causes GPU thread divergence but is bounded.
- **Memory access**: Each streamline step performs 6 bilinear interpolations on 2D
  regular grids. Access pattern is irregular (trajectory-following), but the grid is
  regular so interpolation is simple (no unstructured mesh lookups).
- **Branching**: Latitude/longitude wrapping has 4 conditional branches per step.
  Grid points below precipitation threshold are skipped entirely. Both are manageable
  on GPU.
- **Arithmetic intensity**: Moderate. Each step is ~30 FLOPs plus 6 memory lookups.
  Likely memory-bandwidth-bound on GPU.
- **Data types**: All Float64 in Matlab. Float32 is likely sufficient for the physics
  (δ values are O(10), τ is O(1–10)) but must be validated.

### Dependencies

- **Pure Matlab** — no toolbox dependencies.
- **`inpaint_nans`** (John D'Errico, FileExchange): Used once at setup to fill NaN
  gaps in `delta_e` and `Tcond`. Not performance-critical.
- **CSV file**: Small lookup table of condensation temperature vs. equilibrium
  fractionation factor (from Siler et al., 2021). Read once at setup.

---

## Project Structure

```
DeltapJulia/
├── Project.toml
├── src/
│   ├── DeltapModel.jl          # module definition
│   ├── types.jl                # shared structs and type aliases
│   ├── interpolation.jl        # bilinear interpolation on regular grid
│   ├── inpaint.jl              # NaN infilling (Laplacian iterative)
│   ├── setup.jl                # grid extension, interpolant precomputation
│   ├── streamline_cpu.jl       # Stage 1–2: CPU streamline tracer
│   ├── streamline_gpu.jl       # Stage 3: KernelAbstractions GPU kernel
│   └── postprocess.jl          # per-point reduction (shared logic)
├── test/
│   ├── runtests.jl
│   ├── test_interpolation.jl
│   ├── test_reference.jl       # comparison against Matlab reference data
│   └── fixtures/               # reference .mat or .nc files
├── bench/
│   ├── bench_cpu.jl            # Stage 2 benchmarks
│   ├── bench_gpu.jl            # Stage 3 benchmarks
│   └── compare.jl              # tabular comparison across stages
└── scripts/
    ├── generate_reference.m    # Matlab script for collaborator to run
    └── run_example.jl          # end-to-end example driver
```

---

## Stage 0: Reference Data & Test Harness

### Goal

Establish ground-truth outputs from the original Matlab code for regression testing
at every subsequent stage.

### Approach — collaborator generates reference data

Since the developer does not have a Matlab license, the collaborator should run the
included `generate_reference.m` script (see below) to produce `.mat` files containing
all inputs and outputs for 2–3 test cases at different resolutions.

**Fallback**: If the collaborator is unavailable, attempt to run in GNU Octave (v7+).
The main compatibility risk is `griddedInterpolant` behavior with transposed grid
arguments. Test with a small grid first. Do NOT convert to Python as an intermediate
step — this introduces a second translation that could mask bugs.

### Reference generation script (for collaborator)

Create `scripts/generate_reference.m`:

```matlab
function generate_reference()
% Generate reference input/output files for Julia port validation.
% Run this in Matlab. Produces .mat files in v7 format (HDF5-compatible).

    % --- Case 1: Small diagnostic (e.g., 30x60) ---
    % Load or construct small test inputs here.
    % [Collaborator fills in actual data loading]
    
    tic;
    [deltap_bar, tau_bar] = get_deltap_fast_optimized(E, P, UQ, VQ, Tcond, ...
        LAT, LON, LAT2, LON2, delta_e, csv_file, Plim, dmax);
    elapsed = toc;
    
    save('reference_small.mat', '-v7', ...
        'E', 'P', 'UQ', 'VQ', 'Tcond', 'LAT', 'LON', 'LAT2', 'LON2', ...
        'delta_e', 'Plim', 'dmax', 'deltap_bar', 'tau_bar', 'elapsed');
    
    fprintf('Small case: %.2f seconds\n', elapsed);
    
    % --- Case 2: Medium (e.g., 180x360) ---
    % ... repeat pattern ...
    
    % --- Case 3: Production resolution ---
    % ... repeat pattern ...
end
```

**Important**: The collaborator should also save the CSV file contents (alpha_eq table)
inside the .mat file or as a separate fixture file.

### Test harness

Create `test/test_reference.jl`:

```julia
using Test, MAT

@testset "Reference comparison" begin
    for case_file in readdir("test/fixtures"; join=true)
        endswith(case_file, ".mat") || continue
        
        data = matread(case_file)
        
        # Run Julia implementation
        deltap_bar, tau_bar = get_deltap_fast_optimized(
            data["E"], data["P"], data["UQ"], data["VQ"], data["Tcond"],
            data["LAT"], data["LON"], data["LAT2"], data["LON2"],
            data["delta_e"], alpha_eq_table, data["Plim"], data["dmax"]
        )
        
        ref_dp = data["deltap_bar"]
        ref_tau = data["tau_bar"]
        
        # Compare where both are non-NaN
        valid = .!isnan.(ref_dp) .& .!isnan.(deltap_bar)
        
        @test all(isnan.(deltap_bar[isnan.(ref_dp)]))  # NaN pattern matches
        @test maximum(abs.(deltap_bar[valid] .- ref_dp[valid])) < 1e-4
        @test maximum(abs.(tau_bar[valid] .- ref_tau[valid])) < 1e-4
    end
end
```

**Tolerance notes**:
- Stage 1 (Float64 CPU): Should match to ~1e-10 or better.
- Stage 2 (optimized CPU): Same precision, just faster.
- Stage 3 (Float32 GPU): Relax to ~1e-3 for deltap_bar, ~1e-4 for tau_bar.
  Validate that the relaxation is acceptable for the science.

---

## Stage 1: Naive Julia CPU Port (Correctness Gate)

### Goal

Line-for-line translation from Matlab. No optimization. Prove correctness against
Stage 0 reference data.

### Translation map

| Matlab construct | Julia equivalent | Notes |
|------------------|------------------|-------|
| `griddedInterpolant(LAT', LON', F', 'linear')` | Manual bilinear interpolation function | See `interpolation.jl` spec below |
| `inpaint_nans(X)` | Iterative Laplacian infill | See `inpaint.jl` spec below |
| `cat(1, A(end,:), A, A(1,:))` | `vcat(A[end:end,:], A, A[1:1,:])` | Cyclic extension |
| `load(csv_file)` | `readdlm(csv_file, ',')` or `CSV.read` | Small lookup table |
| `interp1(..., 'linear', 'extrap')` | Manual 1D linear interpolation with extrapolation | ~10 lines |
| `cumsum(x)` | `cumsum(x)` | Direct |
| `nan(A, B)` | `fill(NaN, A, B)` | Direct |
| `isnan(x(:))` | `any(isnan, x)` | Vectorized check |
| `zeros(N, 1)` | `zeros(N)` | Column vectors |
| `cosd(x)` | `cosd(x)` | In Julia base |
| `sqrt(a.^2 + b.^2)` | `hypot(a, b)` | More numerically stable |

### `interpolation.jl` — Bilinear interpolation on regular grid

```julia
"""
    bilinear_interp(field, lat_grid, lon_grid, lat, lon)

Bilinear interpolation on a regular 2D grid. Assumes:
- `lat_grid` is a vector of increasing latitudes (rows of the grid)
- `lon_grid` is a vector of increasing longitudes (columns of the grid)
- `field` is size (length(lat_grid), length(lon_grid))
- Caller handles cyclic longitude wrapping before calling.

Returns interpolated value at (lat, lon).
"""
function bilinear_interp(field::AbstractMatrix, lat_grid::AbstractVector,
                         lon_grid::AbstractVector, lat::Real, lon::Real)
    # Find bracketing indices
    # lat_grid and lon_grid are assumed sorted ascending
    # Use searchsortedfirst for O(log n) lookup
    
    i = clamp(searchsortedfirst(lat_grid, lat) - 1, 1, length(lat_grid) - 1)
    j = clamp(searchsortedfirst(lon_grid, lon) - 1, 1, length(lon_grid) - 1)
    
    # Fractional positions
    t_lat = (lat - lat_grid[i]) / (lat_grid[i+1] - lat_grid[i])
    t_lon = (lon - lon_grid[j]) / (lon_grid[j+1] - lon_grid[j])
    
    # Bilinear blend
    return (1 - t_lat) * (1 - t_lon) * field[i, j] +
           (1 - t_lat) *      t_lon  * field[i, j+1] +
                t_lat  * (1 - t_lon) * field[i+1, j] +
                t_lat  *      t_lon  * field[i+1, j+1]
end
```

**Critical note on array orientation**: Matlab's `griddedInterpolant(LAT', LON', F')`
transposes the grids. In the Matlab code, `LAT` and `LON` are 2D meshgrids where
rows vary in longitude and columns vary in latitude (Matlab convention). The
transpose makes columns vary in longitude. In Julia, work in native column-major:
store fields as `(n_lat, n_lon)` and extract 1D coordinate vectors from the grid.
**Get this mapping right before anything else — it's the #1 source of port bugs.**

### `inpaint.jl` — NaN infilling

```julia
"""
    inpaint_nans!(A::Matrix{T}; max_iter=1000, tol=1e-6) where T

Simple iterative Laplacian infill of NaN values in a 2D matrix.
Replaces NaN entries with the average of their non-NaN neighbors,
iterating until convergence. Modifies A in-place.
"""
function inpaint_nans!(A::Matrix{T}; max_iter=1000, tol=1e-6) where T
    nan_mask = isnan.(A)
    any(nan_mask) || return A
    
    # Initialize NaN locations with mean of non-NaN values
    μ = mean(filter(!isnan, A))
    A[nan_mask] .= μ
    
    for iter in 1:max_iter
        max_change = zero(T)
        for j in axes(A, 2), i in axes(A, 1)
            nan_mask[i, j] || continue
            n = 0; s = zero(T)
            if i > 1;            n += 1; s += A[i-1, j]; end
            if i < size(A, 1);   n += 1; s += A[i+1, j]; end
            if j > 1;            n += 1; s += A[i, j-1]; end
            if j < size(A, 2);   n += 1; s += A[i, j+1]; end
            new_val = s / n
            max_change = max(max_change, abs(new_val - A[i, j]))
            A[i, j] = new_val
        end
        max_change < tol && break
    end
    return A
end
```

### Main function skeleton (Stage 1)

The function signature should closely mirror the Matlab version:

```julia
function get_deltap_cpu(E, P, UQ, VQ, Tcond, LAT, LON, LAT2, LON2,
                        delta_e, alpha_eq, Plim, dmax)
    # --- Setup (mirrors Matlab lines 1–55) ---
    dx = 15.0          # km
    Dx = dx / 111.0    # deg latitude increment
    Nmax = ceil(Int, dmax / dx)
    
    # Clamp delta_e
    delta_e = clamp.(delta_e, -200.0, 200.0)
    
    # Inpaint NaNs (mutating copies)
    delta_e = copy(delta_e); inpaint_nans!(delta_e)
    Tcond = copy(Tcond); inpaint_nans!(Tcond)
    
    # Cyclic grid extension (vcat wrapping)
    # ... [see translation map above]
    
    # Extract 1D coordinate vectors for interpolation
    lat_vec = LAT[:, 1]   # assuming LAT varies along dim 1
    lon_vec = LON[1, :]   # assuming LON varies along dim 2
    # *** VERIFY THIS ORIENTATION AGAINST MATLAB ***
    
    # --- Pre-interpolate P2, Plim2 at output grid ---
    A, B = size(LAT2)
    deltap_bar = fill(NaN, A, B)
    tau_bar    = fill(NaN, A, B)
    
    # --- Main loop ---
    for j in 1:B, i in 1:A      # column-major iteration order
        lat2_ij = LAT2[i, j]
        lon2_ij = LON2[i, j]
        
        P2_ij = bilinear_interp(P, lat_vec, lon_vec, lat2_ij, lon2_ij)
        Plim2_ij = bilinear_interp(Plim, lat_vec, lon_vec, lat2_ij, lon2_ij)
        P2_ij < Plim2_ij && continue
        
        # Streamline tracing (while loop)
        # ... [direct translation of lines 73–115]
        
        # Post-processing (cumsum, alpha lookup, weighting)
        # ... [direct translation of lines 120–155]
    end
    
    return deltap_bar, tau_bar
end
```

### Acceptance criteria

- `test/test_reference.jl` passes at `rtol=1e-8` for all reference cases.
- No GPU code, no threading, no optimization.

---

## Stage 2: Optimized Julia CPU (Threaded)

### Goal

Maximize CPU performance for a fair Matlab-vs-Julia comparison. This is the
CPU↔CPU benchmark stage.

### Optimizations (in priority order)

#### 2a. Thread the outer loop

```julia
using Base.Threads

# Pre-allocate per-thread work arrays
work = [
    (tau=zeros(Nmax), E0=zeros(Nmax), Tcond0=zeros(Nmax),
     Fmag0=zeros(Nmax), P0=zeros(Nmax), delta_e0=zeros(Nmax))
    for _ in 1:nthreads()
]

@threads for idx in 1:(A * B)
    i, j = divrem(idx - 1, A) .+ (1, 1)  # column-major linear index
    tid = threadid()
    w = work[tid]
    # ... use w.tau, w.E0, etc. instead of allocating inside loop
end
```

#### 2b. Precompute 1D coordinate vectors and grid spacing

Instead of calling `searchsortedfirst` per interpolation, precompute:
- `inv_dlat = 1.0 / (lat_vec[2] - lat_vec[1])` (regular grid!)
- `inv_dlon = 1.0 / (lon_vec[2] - lon_vec[1])`

Then index computation is a single multiply + floor:
```julia
i = clamp(floor(Int, (lat - lat_vec[1]) * inv_dlat) + 1, 1, n_lat - 1)
```

This eliminates the binary search entirely — significant for millions of interpolation
calls.

#### 2c. Precompute alpha_eq as a fast 1D lookup

The α_eq table is small. Build a regularly-spaced lookup vector at setup:
```julia
Tmin, Tmax = extrema(alpha_eq[:, 1])
dT = 0.1  # fine enough resolution
T_lookup = range(Tmin, Tmax; step=dT)
alpha_lookup = [interp1d(alpha_eq, T) for T in T_lookup]
```
Then at runtime, the lookup is `floor` + lerp with no search.

#### 2d. Column-major loop ordering

Ensure the outer loop iterates `for j in 1:B, i in 1:A` (Julia column-major).
This matters for the output array writes and for any cache-line effects on input
array reads, though the streamline access pattern is inherently irregular.

#### 2e. Optional: Float32

If Stage 1 passes at Float32 tolerance (~1e-3 for deltap_bar), switch to Float32
here. This halves memory traffic and sets up the GPU stage. Test separately.

### Benchmarking

```julia
# bench/bench_cpu.jl
using BenchmarkTools

# Warmup (JIT compilation)
get_deltap_cpu(small_inputs...)

# Benchmark
t_julia = @belapsed get_deltap_cpu($inputs...)

println("Julia CPU ($(nthreads()) threads): $(t_julia) s")
println("Matlab CPU (from reference):       $(matlab_elapsed) s")
println("Speedup: $(matlab_elapsed / t_julia)×")
```

### Acceptance criteria

- Reference test still passes (identical outputs to Stage 1 at Float64, or within
  Float32 tolerance if using Float32).
- Wall-clock benchmarks captured at ≥2 resolutions.
- Speedup vs Matlab documented.

---

## Stage 3: GPU via KernelAbstractions.jl

### Goal

Port the embarrassingly parallel outer loop to GPU with a single kernel that runs
on both NVIDIA (CUDA) and Apple Silicon (Metal).

### Package dependencies

```toml
[deps]
KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"       # NVIDIA backend
Metal = "dde4c033-4e86-420c-a63e-0dd931031962"       # Apple Silicon backend
Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"       # struct adaptation
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182" # fixed-size arrays in kernel
```

### Backend selection

```julia
function select_backend()
    if CUDA.functional()
        return CUDABackend(), CuArray
    elseif Metal.functional()
        return MetalBackend(), MtlArray
    else
        return CPU(), Array   # fallback to KA's CPU backend
    end
end
```

### Data layout for GPU

Define a struct holding all the read-only input fields on the device:

```julia
struct DeltapFields{T, M <: AbstractMatrix{T}, V <: AbstractVector{T}}
    P::M
    UQ::M
    VQ::M
    E::M
    Tcond::M
    delta_e::M
    Plim::M
    lat_vec::V
    lon_vec::V
    inv_dlat::T
    inv_dlon::T
    lat0::T          # lat_vec[1]
    lon0::T          # lon_vec[1]
    n_lat::Int32
    n_lon::Int32
    alpha_T::V        # pre-sampled alpha_eq temperatures
    alpha_val::V      # pre-sampled alpha_eq values
    alpha_Tmin::T
    alpha_inv_dT::T
end

Adapt.@adapt_structure DeltapFields
```

### GPU kernel

```julia
using KernelAbstractions, StaticArrays

@kernel function deltap_kernel!(deltap_bar, tau_bar, @Const(LAT2), @Const(LON2),
                                 @Const(fields), dx, Dx, Nmax_val)
    idx = @index(Global, Linear)
    A = size(LAT2, 1)
    # Column-major (i, j) from linear index
    j = (idx - 1) ÷ A + 1
    i = (idx - 1) % A + 1
    
    lat2 = LAT2[i, j]
    lon2 = LON2[i, j]
    
    # Quick rejection
    P2_val = bilinear_interp_gpu(fields.P, fields, lat2, lon2)
    Plim2_val = bilinear_interp_gpu(fields.Plim, fields, lat2, lon2)
    if P2_val < Plim2_val
        return
    end
    
    # --- Streamline tracing ---
    # Use MVector for register-resident per-thread arrays.
    # IMPORTANT: Nmax_val must be a compile-time constant for MVector.
    # If Nmax varies, use a fixed upper bound (e.g., 512) and track
    # actual length with a counter.
    
    NMAX = 512  # compile-time upper bound; adjust based on typical dmax
    
    tau   = @MVector zeros(Float32, NMAX)
    E0    = @MVector zeros(Float32, NMAX)
    P0    = @MVector zeros(Float32, NMAX)
    Tc0   = @MVector zeros(Float32, NMAX)
    Fm0   = @MVector zeros(Float32, NMAX)
    de0   = @MVector zeros(Float32, NMAX)
    
    lat_cur = lat2
    lon_cur = lon2
    
    vq_cur = bilinear_interp_gpu(fields.VQ, fields, lat_cur, lon_cur)
    uq_cur = bilinear_interp_gpu(fields.UQ, fields, lat_cur, lon_cur)
    E0[1]  = bilinear_interp_gpu(fields.E, fields, lat_cur, lon_cur)
    Tc0[1] = bilinear_interp_gpu(fields.Tcond, fields, lat_cur, lon_cur)
    Fm0[1] = sqrt(uq_cur^2 + vq_cur^2)
    P0[1]  = P2_val
    de0[1] = bilinear_interp_gpu(fields.delta_e, fields, lat_cur, lon_cur)
    
    jj = 1
    while tau[jj] < 10f0 && jj < min(Nmax_val, NMAX) - 1
        inv_fm = 1f0 / Fm0[jj]
        dtheta = -vq_cur * inv_fm * Float32(Dx)
        dphi   = -uq_cur * inv_fm * Float32(Dx) / cosd(lat_cur + dtheta / 2f0)
        
        lat_new = lat_cur + dtheta
        lon_new = lon_cur + dphi
        
        # Wrapping (same logic as Matlab)
        if lat_new < -90f0
            lat_new = 180f0 + lat_new; lon_new -= 180f0
        elseif lat_new > 90f0
            lat_new = 180f0 - lat_new; lon_new -= 180f0
        end
        if lon_new < 0f0
            lon_new += 360f0
        elseif lon_new > 360f0
            lon_new -= 360f0
        end
        
        vq_cur = bilinear_interp_gpu(fields.VQ, fields, lat_new, lon_new)
        uq_cur = bilinear_interp_gpu(fields.UQ, fields, lat_new, lon_new)
        jj += 1
        E0[jj]  = bilinear_interp_gpu(fields.E, fields, lat_new, lon_new)
        Tc0[jj] = bilinear_interp_gpu(fields.Tcond, fields, lat_new, lon_new)
        Fm0[jj] = sqrt(uq_cur^2 + vq_cur^2)
        P0[jj]  = bilinear_interp_gpu(fields.P, fields, lat_new, lon_new)
        de0[jj] = bilinear_interp_gpu(fields.delta_e, fields, lat_new, lon_new)
        
        mu_val = P0[jj-1] / Fm0[jj-1]
        tau[jj] = tau[jj-1] + mu_val * Float32(dx) * 1000f0
        
        lat_cur = lat_new
        lon_cur = lon_new
    end
    
    npts = jj  # actual streamline length
    
    # --- Post-processing (serial, per-thread) ---
    # Recompute tau via cumsum (matches Matlab post-loop recomputation)
    # ... [see post-processing section below]
    
    # Write results
    # deltap_bar[i, j] = ...
    # tau_bar[i, j] = ...
end
```

### GPU bilinear interpolation

```julia
@inline function bilinear_interp_gpu(field, fields, lat, lon)
    i = clamp(floor(Int32, (lat - fields.lat0) * fields.inv_dlat) + Int32(1),
              Int32(1), fields.n_lat - Int32(1))
    j = clamp(floor(Int32, (lon - fields.lon0) * fields.inv_dlon) + Int32(1),
              Int32(1), fields.n_lon - Int32(1))
    
    t_lat = (lat - fields.lat_vec[i]) * fields.inv_dlat
    t_lon = (lon - fields.lon_vec[j]) * fields.inv_dlon
    
    return (1 - t_lat) * (1 - t_lon) * field[i, j] +
           (1 - t_lat) *      t_lon  * field[i, j+1] +
                t_lat  * (1 - t_lon) * field[i+1, j] +
                t_lat  *      t_lon  * field[i+1, j+1]
end
```

### GPU post-processing

The Matlab code's post-processing uses `cumsum`, `interp1`, `exp`, `flip`, and
conditional indexing. On GPU, all of this must be done per-thread with serial loops
(no Julia stdlib calls available inside `@kernel`):

```julia
# Inside the kernel, after streamline tracing:

# 1. Recompute tau = cumsum(P0 ./ Fm0) .* dx * 1000
#    (serial cumsum over npts elements)
running = 0f0
for k in 1:npts
    running += P0[k] / Fm0[k]
    tau[k] = running * Float32(dx) * 1000f0
end

# 2. Alpha_eq lookup (1D interpolation from pre-sampled table)
#    For each k: alphac0[k] = lookup(Tc0[k] - 273.15)
#    Use the pre-sampled regular table in fields.alpha_T / fields.alpha_val

# 3. Valid mask: tau > 0 && -120 < de0 < 120; always include k=1

# 4. Exponential decay weighting (reverse cumsum of mu * dx * 1000)
#    This requires a backward pass — compute from npts down to 1

# 5. Weighted mean for alpha_bar, then delta_p, then final weighted average
```

### Launch configuration

```julia
function run_gpu(E, P, UQ, VQ, Tcond, LAT, LON, LAT2, LON2,
                 delta_e, alpha_eq, Plim, dmax)
    backend, ArrayType = select_backend()
    
    # Setup (same as CPU: clamp, inpaint, cyclic extension)
    # ... then move to device:
    fields = DeltapFields(
        ArrayType(P_ext), ArrayType(UQ_ext), ...
    )
    
    d_LAT2 = ArrayType(LAT2)
    d_LON2 = ArrayType(LON2)
    A, B = size(LAT2)
    d_deltap = ArrayType(fill(NaN32, A, B))
    d_tau    = ArrayType(fill(NaN32, A, B))
    
    Nmax = ceil(Int, dmax / 15.0)
    
    kernel! = deltap_kernel!(backend, 256)  # workgroup size 256
    kernel!(d_deltap, d_tau, d_LAT2, d_LON2, fields,
            Float32(15.0), Float32(15.0/111.0), Int32(Nmax);
            ndrange=A * B)
    KernelAbstractions.synchronize(backend)
    
    return Array(d_deltap), Array(d_tau)
end
```

### Acceptance criteria

- Reference test passes at Float32 tolerance (rtol ≈ 1e-3).
- Runs correctly on both NVIDIA GPU and Apple M-series GPU.
- Speedup vs Stage 2 CPU documented at ≥2 resolutions.

---

## Stage 4: Profiling & Tuning

### NVIDIA (Nsight Compute / CUDA.@profile)

Key metrics to check:
- **Occupancy**: If < 50%, likely register pressure from large `MVector`s. Consider
  reducing NMAX or splitting into smaller kernels.
- **Memory throughput**: Compare achieved bandwidth to peak. If < 50% of peak,
  consider texture memory for read-only fields (CUDA-specific optimization, would
  need a separate CUDA kernel path).
- **Warp divergence**: Profile the while-loop termination. If some grid regions
  consistently use few steps, consider sorting grid points by expected path length.
- **Register spill to local memory**: Check if `MVector` arrays spill. If so,
  reduce per-thread array sizes or use a two-pass approach.

### Apple Silicon (Xcode Instruments GPU profiler)

Key metrics:
- **ALU utilization**: M-series GPUs have high ALU:bandwidth ratios.
- **Threadgroup occupancy**: Metal has different limits than CUDA (max 1024
  threads per threadgroup, but register pressure differs).
- **Bandwidth**: Check if the irregular memory access pattern (trajectory-following
  interpolation) is cache-friendly. M-series unified memory eliminates PCIe transfer
  overhead, which is an advantage for this workload.

### Tuning knobs

| Parameter | Starting value | Tuning range | Effect |
|-----------|---------------|-------------|--------|
| Workgroup size | 256 | 64–1024 | Occupancy vs register pressure |
| NMAX (compile-time) | 512 | 128–1024 | Register pressure; must cover dmax |
| Float type | Float32 | Float32/Float64 | 2× bandwidth + throughput on consumer GPU |
| Grid-point sorting | Off | On/off | Reduces warp divergence if path lengths vary |

### Benchmarking matrix

```
bench/compare.jl should produce a table like:

| Resolution | Matlab CPU | Julia CPU 1T | Julia CPU 8T | Julia GPU (NVIDIA) | Julia GPU (Metal) |
|------------|-----------|-------------|-------------|-------------------|-------------------|
| 30×60      |    X.Xs   |     X.Xs    |    X.Xs     |       X.Xs        |       X.Xs        |
| 180×360    |    X.Xs   |     X.Xs    |    X.Xs     |       X.Xs        |       X.Xs        |
| 720×1440   |    X.Xs   |     X.Xs    |    X.Xs     |       X.Xs        |       X.Xs        |
```

---

## Known Risks and Mitigations

### 1. Array orientation bug (HIGH probability)

The Matlab code uses `griddedInterpolant(LAT', LON', F')` with transposed inputs.
Matlab meshgrids and Julia array conventions differ. **Mitigation**: Build a tiny
3×4 synthetic test case with known analytic values and verify interpolation matches
Matlab output before proceeding to the full port.

### 2. Metal.jl maturity

Metal.jl via KernelAbstractions is newer than CUDA.jl and may have edge cases
(e.g., `StaticArrays` MVector support, `cosd` availability on Metal). **Mitigation**:
Start GPU work on NVIDIA where the tooling is mature, then port to Metal. If Metal
lacks `cosd`, use `cos(lat * π / 180f0)`.

### 3. MVector register spill

If `NMAX` is large (e.g., >300), six `MVector{NMAX, Float32}` per thread may
exceed register limits and spill to slow local memory. **Mitigation**: Profile early.
If spilling, consider a two-pass approach: (1) trace streamlines and store intermediate
results in global memory, (2) run post-processing as a separate kernel.

### 4. Thread divergence from variable-length streamlines

Grid points over ocean vs. continental interiors may have very different streamline
lengths. **Mitigation**: Acceptable for a first pass. If profiling shows >30%
warp idle time, implement grid-point sorting by estimated path length (e.g., based
on local P/IVT ratio) to group similar-length streamlines into the same warps.

### 5. `inpaint_nans` equivalence

The Matlab `inpaint_nans` uses a sparse matrix solver (spring metaphor). The
iterative Laplacian proposed here is simpler and may give slightly different results
in NaN-dense regions. **Mitigation**: Compare infilled fields directly against Matlab
output. If divergence is significant, implement a closer approximation (e.g., using
`SparseArrays` + `\` solve) or have the collaborator save pre-infilled fields.

---

## Appendix A: Original Matlab Source

See uploaded file `get_deltap_fast_optimized.m` — the complete source is the
authoritative reference for all translation decisions.

## Appendix B: Key Julia Packages & Versions

| Package | Purpose | Min version |
|---------|---------|-------------|
| `KernelAbstractions.jl` | Backend-agnostic GPU kernels | ≥ 0.9 |
| `CUDA.jl` | NVIDIA backend | ≥ 5.0 |
| `Metal.jl` | Apple Silicon backend | ≥ 1.0 |
| `StaticArrays.jl` | Fixed-size arrays in GPU kernels | ≥ 1.9 |
| `MAT.jl` | Read .mat reference files | ≥ 0.10 |
| `BenchmarkTools.jl` | Timing | ≥ 1.3 |
| `DelimitedFiles` (stdlib) | Read CSV alpha_eq table | stdlib |
| `Test` (stdlib) | Test harness | stdlib |
