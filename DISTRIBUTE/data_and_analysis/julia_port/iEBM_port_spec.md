# iEBM Port Specification — Claude Code Task Sheet

## Project Summary

Port the isotope-enabled Energy Balance / Attenuation Model (iEBM) from Python (~2,700 LOC) to a compiled language with GPU acceleration. The model computes τ̄ (tau-bar), the evaporation-weighted mean isotopic attenuation parameter, by marching moisture transport streamlines upstream through gridded climatological fields. Current Python runtime is on the order of **days** for full decomposition sweeps and will scale to **global domains** (1–10 GB NetCDF inputs).

---

## Language Recommendation: Julia

### Rationale

| Factor | Julia | C++ | Winner |
|--------|-------|-----|--------|
| GPU portability (Metal + CUDA) | `Metal.jl` + `CUDA.jl` via `KernelAbstractions.jl` — write once, run on both | Kokkos or SYCL — mature but heavy build systems, no first-class Metal | **Julia** |
| Scientist maintainability | Math-like syntax, REPL-driven development, no build system | CMake/Kokkos boilerplate, longer dev cycles | **Julia** |
| NetCDF I/O | `NCDatasets.jl` (direct, no xarray-like layer needed) | `netcdf-cxx4` — works but verbose | **Julia** |
| Prototyping speed | Interactive; test kernels in minutes | Compile-test cycles much slower | **Julia** |
| Performance ceiling | Within ~5–15% of C++ for this workload pattern | Marginally higher raw throughput | Tie |
| EBM extensibility | Easy to add new physics modules, hot-reload | Requires recompilation | **Julia** |
| Ecosystem fit | Strong in climate/geoscience community | Strong in HPC but heavier | **Julia** |

**The one risk with Julia:** `Metal.jl` is newer and less battle-tested than `CUDA.jl`. Mitigation: the M2 Ultra's unified memory architecture means CPU-side Julia is already fast on Apple Silicon (no discrete GPU memory transfer overhead), so Metal acceleration is a bonus rather than a requirement. The code will be structured so the Metal backend can mature while NVIDIA HPC runs remain the primary GPU target.

---

## Architecture Overview

```
iEBM.jl/
├── Project.toml                  # Dependencies
├── src/
│   ├── iEBM.jl                   # Module root, exports
│   ├── Config.jl                 # Run configuration (replaces initialize.py dict)
│   ├── IO.jl                     # NetCDF read/write, coordinate setup
│   ├── Hydroclim.jl              # Orographic partitioning, field prep
│   ├── EBM.jl                    # Energy balance model (stub + hooks for reactivation)
│   ├── Interpolation.jl          # RegularGridInterpolator equivalent (GPU-compatible)
│   ├── Streamlines.jl            # Core streamline marching kernel (HOT PATH)
│   ├── TauBar.jl                 # τ̄ integration, E-weighting, land fraction
│   ├── Decomposition.jl          # Spatial, climatological, E/L/path decomposition
│   └── GPUKernels.jl             # KernelAbstractions.jl wrappers
├── test/
│   ├── runtests.jl
│   ├── unit/                     # Per-function tests against analytic/hand-calculated values
│   │   ├── test_interpolation.jl
│   │   ├── test_streamline_march.jl
│   │   ├── test_taubar_integration.jl
│   │   ├── test_terrestrial_efrac.jl
│   │   ├── test_decomposition.jl
│   │   ├── test_hydroclim.jl
│   │   └── test_config.jl
│   ├── integration/              # Full-pipeline runs on small inputs
│   │   ├── test_houston_bbox.jl
│   │   ├── test_coord_list.jl
│   │   └── test_decomposition_run.jl
│   ├── regression/               # Versioned golden-reference framework
│   │   ├── RegressionRunner.jl
│   │   ├── cases/                # Auto-discovered, one dir per case
│   │   │   ├── houston_bbox/
│   │   │   ├── coord_list_basic/
│   │   │   └── ...               # New cases added as model grows
│   │   └── reports/
│   ├── gpu/
│   │   ├── test_cpu_gpu_agreement.jl
│   │   └── test_f32_precision.jl
│   └── fixtures/                 # Small NetCDF extracts for CI
├── benchmark/
│   ├── bench_streamline.jl       # Single-point profiling
│   └── bench_fullrun.jl          # End-to-end scaling
└── scripts/
    ├── run.jl                    # Entry point (replaces RUN.py)
    └── validate_vs_python.jl     # Cross-check tool
```

---

## Module-by-Module Conversion Plan

### 1. `Config.jl` — Run Configuration

**Source:** `initialize.py` (161 lines)

Replace the Python dictionary with a Julia struct hierarchy using `Base.@kwdef` for keyword construction and defaults. Group parameters into nested structs for clarity.

```
RunConfig
├── GridConfig         (deg_per_lat, deg_per_lon)
├── VariableNames      (precip, spechum, evap, UQ, VQ, ...)
├── TopoConfig         (orog_partition, slope_threshold, elev_threshold)
├── ForcingConfig      (insolation, albedo, arbitrary forcing params)
├── EBMConfig          (EFE thresholds, EFPM params, heaviside limits)
├── IsotopeConfig      (tau-bar bbox/coord_list, dx, taumax, dmax, n_samples)
├── DecompConfig       (local_v_regional flags, E/L/s flags, thresholds)
├── IOConfig           (run_name, run_path, input_dir, filenames)
└── Constants          (Lv, g)
```

**Key design choice:** Make the config immutable (`struct`, not `mutable struct`) and pass it as a const-ref throughout the pipeline. This enables the compiler to propagate constants into hot loops.

**Validation:** Add a `validate(config::RunConfig)` function that checks parameter consistency at startup (e.g., bbox ranges, resolution compatibility).


### 2. `IO.jl` — Data I/O

**Source:** `Run.inputfiles()` in `attenuationMod_fxns.py` (lines 1441–1498), plus all `xr.open_dataset` / `to_netcdf` calls.

Replace xarray with `NCDatasets.jl` for reading and writing. Internally, store fields as plain Julia arrays (`Matrix{Float64}` or `Array{Float64,3}`) rather than xarray-like labeled structures. Coordinate metadata lives in the `RunConfig`.

**Functions:**
- `load_climatology(config) → ClimatologyData` — reads the clim NetCDF, returns a struct of arrays
- `load_forcing(config) → Union{ForcingData, Nothing}` — reads forcing or returns nothing
- `load_coordinates(config) → Union{DataFrame, Nothing}` — reads coord list CSV if needed
- `save_results(outdata, streamlines, config)` — writes NetCDF + CSV output

**Data struct pattern:**
```julia
struct ClimatologyData{T <: AbstractFloat, A <: AbstractMatrix{T}}
    lat::Vector{T}
    lon::Vector{T}
    P::A          # precipitation [lat × lon]  (or [lat × lon × slice])
    E::A          # evaporation
    UQ::A         # zonal moisture flux
    VQ::A         # meridional moisture flux
    Q::A          # specific humidity
    LANDFRAC::A   # land fraction
    # ... etc
end
```

This struct is generic over the array type `A`, so it can hold `Matrix{Float64}` on CPU or `CuMatrix{Float64}` / `MtlMatrix{Float64}` on GPU without code changes.


### 3. `Hydroclim.jl` — Orographic Partitioning

**Source:** `Hydroclim.orog_partition()` in `attenuationMod_fxns.py` (lines 43–166)

This is straightforward array math: gradient computation, masking, zonal mean subtraction. Port as broadcasted Julia operations. Not a performance bottleneck.

**Functions:**
- `compute_slope(elevation, lat, lon) → slope_matrix`
- `build_orog_mask(slope, elevation, config) → BitMatrix`
- `partition_field!(dataset, field_name, mask, config)` — modifies dataset in-place


### 4. `EBM.jl` — Energy Balance Model (Stub)

**Source:** Commented-out calls to `EBM.divMSE`, `EBM.EFE_calc_heaviside`, `EBM.stitch_efpm_contours`, `EBM.efpm_rank`, `EBM.efpm_contours_to_mask`

**Status:** Code not provided but referenced in the Python. Create module with:
- Stub functions that pass through data unchanged (current behavior)
- Documented interface contracts for future implementation
- A `config.ebm_enabled::Bool` flag that activates/deactivates the EBM path

**Interface to preserve:**
```julia
# These are no-ops when ebm_enabled=false, but define the contract:
function solve_ebm!(data::ClimatologyData, forcing, config::EBMConfig) end
function compute_efe!(data, config) end
function compute_efpm!(data, config) end
```


### 5. `Interpolation.jl` — Fast Field Interpolation

**Source:** `tau_interpolatorInitialize()` and all `RegularGridInterpolator` usage in `accessory_fxns.py` (lines 227–253)

This is **critical for GPU portability**. SciPy's `RegularGridInterpolator` has no GPU equivalent. Implement a custom bilinear interpolator on regular lat/lon grids.

**Design:**
```julia
struct RegularGridInterp{T, A <: AbstractMatrix{T}}
    lat::Vector{T}
    lon::Vector{T}    # extended with wraparound padding
    data::A           # [lat × lon] with wraparound columns
    dlat::T           # grid spacing (assumed uniform)
    dlon::T
    lat0::T           # minimum latitude
    lon0::T           # minimum longitude
end

function (interp::RegularGridInterp)(lat, lon)
    # Bilinear interpolation — 4 lookups, GPU-compatible
    # No branching, no allocations
end
```

**Wraparound handling:** Pre-pad the longitude dimension (as the Python code already does with `np.concatenate`) so the interpolator never needs bounds-checking logic.

**GPU note:** This interpolator is called ~1,800× per streamline × thousands of grid cells. It MUST be allocation-free and branchless for GPU kernels. The struct-of-arrays layout above enables coalesced memory access when many streamlines query nearby points.


### 6. `Streamlines.jl` — Core Marching Kernel ⚡ HOT PATH

**Source:** `tau_streamline_point()` in `accessory_fxns.py` (lines 56–176)

This is where ~95% of compute time is spent. The function marches upstream along moisture transport vectors, interpolating 6 fields at each step and accumulating tau via trapezoidal integration.

**CPU version:**
```julia
function march_streamline!(
    result::StreamlineResult,
    lat0::Float64, lon0::Float64,
    Nmax::Int, taumax::Float64, Dx::Float64, dx::Float64,
    Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit
)
    # Pre-allocate result arrays in StreamlineResult struct
    # Sequential march — each step depends on previous position
    @inbounds for step in 1:Nmax-1
        result.tau[step] >= taumax && break
        
        # Direction from transport field
        dtheta = -result.vq[step] / result.Fmag[step] * Dx
        # Pole clamping
        lat1 = clamp(lat0 + dtheta, -89.5, 89.5)
        dtheta = lat1 - lat0
        dphi = -result.uq[step] / result.Fmag[step] * Dx / cos(deg2rad(lat0 + dtheta/2))
        lon1 = mod(lon0 + dphi, 360.0)
        
        # Interpolate all fields at new position
        result.E[step+1] = Efit(lat1, lon1)
        result.P[step+1] = Pfit(lat1, lon1)
        # ... etc
        
        # Trapezoidal tau integration
        result.tau[step+1] = result.tau[step] + 
            (dx * 1000) * (result.P[step]/result.Fmag[step] + result.P[step+1]/result.Fmag[step+1]) / 2
        
        lat0, lon0 = lat1, lon1
    end
end
```

**GPU parallelization strategy:**

The streamline march is **inherently sequential** within a single grid cell — each step's position depends on the previous. However, all grid cells are **embarrassingly parallel**. The GPU strategy is:

1. **One thread per sink point** (lat/lon grid cell where τ̄ is computed)
2. Each thread runs the full sequential march independently
3. All interpolation data (P, E, UQ, VQ, Fmag, LANDFRAC grids) lives in GPU global memory (fits in M2 Ultra's 192 GB unified memory or NVIDIA HBM)
4. Per-thread workspace (tau, E0, P0, etc. arrays of length Nmax≈1800) allocated as thread-local

**Expected parallelism:**
- Regional 4° bbox at 0.25° resolution: 16 × 16 = 256 parallel streamlines
- Global land at 1° resolution: ~15,000 parallel streamlines
- Global land at 0.25° resolution: ~250,000 parallel streamlines

This maps very well to GPU thread blocks, especially at the global scale. The M2 Ultra has 76 GPU cores; NVIDIA A100 has 6,912 CUDA cores. Even at regional scale (256 threads), this saturates the M2; at global scale it fully saturates an A100.

**KernelAbstractions.jl kernel sketch:**
```julia
@kernel function taubar_kernel!(
    tau_bar, tau_bar_Ewtd, moisture_dist, land_frac_stream, land_frac_Esource,
    lat_grid, lon_grid, lat_bool, lon_bool,
    E_grid, P_grid, UQ_grid, VQ_grid, Fmag_grid, lfrac_grid,
    config
)
    idx = @index(Global)
    # Map linear index → (y, x) in the bbox
    # ... run march_streamline! for this point
    # ... compute tau_bar via trapezoidal integration
    # ... store results
end
```


### 7. `TauBar.jl` — τ̄ Integration

**Source:** `taubar_noDecomposition()` (lines 184–552), `dist_to_tau_coords()`, `terrestrial_E_frac()` in `accessory_fxns.py`

This wraps the streamline kernel: sets up the interpolation grids, dispatches over the lat/lon bbox (or coord list), collects results, and regridds back to the original resolution if needed.

**Functions:**
- `compute_taubar!(output, data, config)` — main entry point
- `taubar_at_point(streamline_result) → Float64` — the τ̄ integration formula
- `terrestrial_efrac(streamline_result, config) → (tau_wtd_mean, land_frac, E_terrFrac)`
- `regrid_to_original!(output, coarse_result, original_coords)` — bilinear regridding if bbox resolution differs from input

**Key formula (line 396 / 510 in Python):**
```
τ̄ = ∫ E₀(τ) · (1/μ(τ)) · τ · e^(-τ) dτ  /  ∫ E₀(τ) · (1/μ(τ)) · e^(-τ) dτ
```

Use `trapz`-equivalent in Julia (manual trapezoidal loop or a thin helper — avoid allocating with `diff()`).


### 8. `Decomposition.jl` — Attribution Analysis

**Source:** `taubar_withDecomposition()` (lines 557–1435), `tau_spatial_decompose()`, `dist_to_Dtaubar_frac_LocalEvap()`, `dist_to_Dtaubar_frac_streamline()`, `tau_decompose_ELPath()`, `local_plus_regional_tauStreamline()` in `accessory_fxns.py`

This is the most complex module and the reason runs take days. It runs 3× more streamlines per cell (for E/L/path attribution) and includes binary search loops for distance-to-fraction computations.

**Functions:**
- `decompose_taubar!(output, data, config)` — main decomposition entry point
- `spatial_decompose!(...)` — local vs. regional attribution
- `ELpath_decompose!(...)` — evaporation, length-scale, flow-path attribution
- `find_fraction_distance(...)` — binary search for distance at which fraction of Δτ̄ is explained

**GPU note:** The decomposition multiplies compute by 3–4× but the parallelism structure is identical (one thread per grid cell). The binary search in `dist_to_Dtaubar_frac_*` functions adds variable-length work per thread (warp divergence on GPU), but this is acceptable given the coarse parallelism granularity.


### 9. `GPUKernels.jl` — Backend Abstraction

**Source:** N/A (new)

Provides the `KernelAbstractions.jl` wrappers and backend selection logic.

```julia
function select_backend(config)
    if CUDA.functional()
        return CUDABackend()
    elseif Metal.functional()
        return MetalBackend()
    else
        return CPU()
    end
end

function to_device(data::ClimatologyData, backend)
    # Convert all arrays to device arrays
    # Leverages the generic A parameter in ClimatologyData
end
```

---

## Phased Implementation Roadmap

### Phase 1 — CPU Julia Port + Validation (Target: 2–3 weeks)

**Goal:** Exact numerical reproduction of Python output on CPU.

1. Implement `Config.jl` with all parameters from `initialize.py`
2. Implement `IO.jl` with NCDatasets.jl
3. Implement `Interpolation.jl` — the custom bilinear interpolator
4. Implement `Streamlines.jl` — `march_streamline!` with pre-allocated buffers
5. Implement `TauBar.jl` — `compute_taubar!` with the bbox loop
6. Implement `Hydroclim.jl` — orographic partitioning
7. Implement `EBM.jl` — stubs only
8. Set up the regression framework: create `RegressionRunner.jl`, register the Houston bbox case, generate the first golden reference from the Python code (`generator = "python"`)
9. Write unit tests for each module against analytic/hand-calculated values
10. Run regression suite in `strict` mode — Julia output must match Python golden reference to <1e-10 relative error
11. Add `Decomposition.jl` — full decomposition paths
12. Register decomposition regression case, generate golden reference, validate

**Expected CPU speedup from Python:** 10–50× from eliminating interpreter overhead, enabling SIMD on the interpolation, and pre-allocating all buffers. This alone may bring days → hours.

**Validation criterion:** All regression cases pass in `strict` mode. When a case fails, the auto-generated report identifies the worst-offending grid cells with lat/lon coordinates for targeted debugging.

### Phase 2 — CPU Performance Optimization (Target: 1 week)

1. **Thread parallelism:** `Threads.@threads` over the (lat, lon) bbox loop — instant N× on multi-core (M2 Ultra has 24 cores → ~20× for compute-bound work)
2. **Pre-allocated workspace pool:** One `StreamlineResult` buffer per thread, reused across grid cells
3. **Interpolation optimization:** Precompute lat/lon index offsets, cache grid spacing reciprocals
4. **Profile with `@btime` and `Profile.jl`** — identify any remaining hotspots
5. **SIMD hints:** `@simd` on the trapezoidal integration loop

**Expected total speedup from Python:** 100–500× (Julia compilation × threading × allocation elimination). Days → minutes for regional runs.

### Phase 3 — GPU Acceleration (Target: 2–3 weeks)

1. Add `KernelAbstractions.jl`, `CUDA.jl`, `Metal.jl` dependencies
2. Implement `GPUKernels.jl` — backend selection, data transfer
3. Write the `taubar_kernel!` — one thread per grid cell
4. Handle the variable-length streamline march (thread divergence) — profile and optimize
5. GPU-compatible interpolation (no allocations, no exceptions)
6. Validate GPU output against CPU output (should be bitwise identical for Float64, near-identical for Float32)
7. Benchmark on M2 Ultra (Metal) and NVIDIA GPU (CUDA)
8. Add Float32 option for additional GPU throughput (with validation of acceptable precision loss)

**Expected GPU speedup over Phase 2 CPU:** 5–20× depending on domain size and hardware. Global domains at 0.25° will benefit most.

### Phase 4 — Production Hardening (Target: 1 week)

1. CLI interface (`run.jl` with argument parsing via `Comonicon.jl` or `ArgParse.jl`)
2. Logging and progress reporting (replace Python print-based tracker)
3. Checkpoint/restart for long decomposition runs (serialize intermediate state to NetCDF)
4. CI with GitHub Actions — unit tests, `strict` regression suite, and GPU agreement tests (see CI config in Regression Testing Framework section)
5. Documentation (Documenter.jl)
6. Package registration (if desired)

---

## Key Technical Decisions

### 1. Float64 vs Float32

Default to **Float64** for scientific accuracy. The τ̄ integration involves exponentials of negative values (e^(-τ) where τ can reach 8), and the trapezoidal sums can accumulate rounding errors. Offer Float32 as an opt-in for GPU runs where speed matters more than precision, with a validation pass to quantify the error.

### 2. Memory Layout

Store all 2D fields as **column-major** `Matrix{Float64}` with dimensions `[lat, lon]` (Julia's natural layout). The interpolator will access `[lat, lon]` pairs, so column-major with lat as the fast index gives good spatial locality for meridional streamlines.

For GPU: the same layout works because each thread accesses its own spatial neighborhood. No need for structure-of-arrays transformation.

### 3. Streamline Workspace

Pre-allocate a `StreamlineResult` struct with fixed-size arrays of length `Nmax`:

```julia
struct StreamlineResult{T}
    tau::Vector{T}      # length Nmax
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
    nsteps::Base.RefValue{Int}  # actual number of steps taken
end
```

On GPU, these become thread-local arrays (either stack-allocated if small enough, or carved out of shared/global memory).

### 4. Coordinate Handling

Drop xarray's labeled-array paradigm. Store coordinates as plain vectors in the config. The Python code uses `ds.sel(lat=..., method='nearest')` extensively — replace with `argmin(abs.(lat .- target))` or a pre-built lookup table.

### 5. Streamline Output (DataFrame)

The Python code builds a DataFrame incrementally with `pd.concat` inside the hot loop (very slow). Replace with:
- Pre-allocate output arrays sized to (n_cells × max_streamline_length / coarsener)
- Fill in-place during the loop
- Convert to DataFrame once at the end via `DataFrames.jl`

---

## Dependencies

```toml
[deps]
NCDatasets = "..."         # NetCDF I/O
DataFrames = "..."         # Streamline output
CSV = "..."                # CSV I/O for coord lists
Interpolations = "..."     # Reference impl (custom will replace in hot path)
StaticArrays = "..."       # Small fixed-size buffers
KernelAbstractions = "..."  # GPU abstraction
CUDA = "..."               # NVIDIA backend
Metal = "..."              # Apple Silicon backend

[extras]
BenchmarkTools = "..."     # Profiling
Test = "..."               # Unit tests
Profile = "..."            # CPU profiling
```

---

## Risk Register

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| `Metal.jl` immaturity / bugs | GPU path broken on Mac | Medium | CPU Julia on M2 Ultra is already fast (24 cores); Metal is a bonus, not a blocker |
| Numerical divergence from Python | Validation failure | Low | Use Float64 throughout; compare at each module boundary via regression framework in `strict` mode |
| EBM reactivation requires major refactor | Scope creep | Medium | Clean interface contracts in `EBM.jl` stubs; pass-through design makes activation additive |
| Variable-length streamlines cause GPU warp divergence | GPU underperformance | Medium | Profile; consider padding short streamlines or sorting cells by expected path length |
| Thread-local memory pressure on GPU | Kernel launch failure | Low | `Nmax ≈ 1800 × 11 arrays × 8 bytes ≈ 155 KB` per thread — fits in GPU registers/local mem for moderate occupancy |
| Physics changes silently break output | Science errors go unnoticed | Medium | Regression suite in `science` mode on every PR; golden reference re-generation requires explicit action and produces a diff report |
| Golden references become stale | False confidence | Low | Manifests record generator commit hash; CI warns if golden ref was generated >N commits behind current HEAD |

---

## Regression Testing Framework

Since the model is under active development, the regression framework must evolve alongside the science. It is designed around three principles: (1) golden references are **versioned and rebuildable**, not static artifacts; (2) tests operate at **multiple granularities** so a failing integration test can be traced to a specific function; (3) new test cases can be added with **minimal friction** as the model grows.


### Directory Structure

```
test/
├── runtests.jl                       # Entry point — runs all suites
├── unit/
│   ├── test_interpolation.jl         # Bilinear interpolator vs analytic functions
│   ├── test_streamline_march.jl      # Single-point march vs hand-calculated values
│   ├── test_taubar_integration.jl    # τ̄ formula on synthetic exponential profiles
│   ├── test_terrestrial_efrac.jl     # Land E-fraction on known geometries
│   ├── test_decomposition.jl         # Spatial / E-L-path decomposition identities
│   ├── test_hydroclim.jl             # Orog partition: mask + zonal mean subtraction
│   └── test_config.jl                # Config validation, defaults, type stability
├── integration/
│   ├── test_houston_bbox.jl          # Full run on Houston bbox case
│   ├── test_coord_list.jl            # Coord-list mode end-to-end
│   └── test_decomposition_run.jl     # Full decomposition sweep
├── regression/
│   ├── RegressionRunner.jl           # Harness: generate, compare, update, report
│   ├── cases/                        # One subdirectory per registered test case
│   │   ├── houston_bbox/
│   │   │   ├── case.toml             # Case metadata + tolerance spec
│   │   │   ├── generate_reference.jl # Script to produce golden output (Python or Julia)
│   │   │   ├── golden/               # Archived reference outputs (git-LFS or .gitignore'd)
│   │   │   │   ├── v001_CLIM.nc
│   │   │   │   ├── v001_STREAMLINES.csv
│   │   │   │   └── v001_manifest.toml  # Hash, date, code version, generator
│   │   │   └── compare.jl            # Case-specific comparison logic (if needed)
│   │   ├── coord_list_basic/
│   │   │   ├── case.toml
│   │   │   ├── ...
│   │   └── global_coarse/            # Added later when global runs become feasible
│   │       ├── case.toml
│   │       └── ...
│   └── reports/                      # Auto-generated comparison reports (.md)
├── gpu/
│   ├── test_cpu_gpu_agreement.jl     # GPU output == CPU output (bitwise for F64)
│   └── test_f32_precision.jl         # Float32 vs Float64 error characterization
└── benchmark/
    ├── bench_streamline.jl           # Single-point profiling
    ├── bench_fullrun.jl              # End-to-end scaling
    └── bench_history.toml            # Tracked performance over time
```


### Test Case Registration (`case.toml`)

Each regression case is self-describing. Adding a new case means creating a directory with a `case.toml` and a reference-generation script — nothing else needs to change.

```toml
[case]
name = "houston_bbox"
description = "Houston region, 0.25° bbox, no decomposition, no forcing"
created = "2026-03-19"
science_version = "v0.1"       # tracks which version of the MODEL PHYSICS this reference reflects
status = "active"               # active | retired | developing

[input]
config_overrides = """
    run_name = "houston"
    Tau-bar_bbox_LatRange = [28, 32]
    Tau-bar_bbox_LonRange = [263, 267]
    Tau-bar_bbox_Resolution = [0.25, 0.25]
    SolveIsotopes = true
    Tau-bar_decomp_E_L_s = false
    Tau-bar_decomp_local_v_regional_LocalEvapEffect = false
    Tau-bar_decomp_local_v_regional_UpwindEffect = false
"""
input_data = "test/fixtures/era_mon_fixvars_houston.nc"   # small extract for CI

[reference]
current_version = "v001"
generator = "python"            # "python" or "julia" — which code produced the golden output
generator_commit = "abc1234"    # git commit of the code that generated this reference

[tolerances]
# Per-variable tolerances: rtol = relative, atol = absolute
# "strict" = must match to floating-point precision (for port validation)
# "science" = must match within scientifically acceptable bounds (for ongoing development)
mode = "strict"

[tolerances.strict]
tau_bar.rtol = 1e-10
tau_bar.atol = 1e-14
tau_bar_wtdLandEvap.rtol = 1e-10
moisture_dist_inland.rtol = 1e-8
streamline_frac_land.rtol = 1e-10
Esource_frac_land.rtol = 1e-10

[tolerances.science]
tau_bar.rtol = 0.01             # 1% — for when physics changes are intentional
tau_bar_wtdLandEvap.rtol = 0.01
moisture_dist_inland.rtol = 0.05
streamline_frac_land.rtol = 0.02
Esource_frac_land.rtol = 0.02
```


### Regression Harness (`RegressionRunner.jl`)

The harness supports four operations:

```julia
module RegressionRunner

using NCDatasets, TOML, Dates, SHA

"""
    run_regression(case_dir; mode=:strict, update=false)

Run a single regression case.
- `mode=:strict` — port validation (tight tolerances)
- `mode=:science` — physics evolution (relaxed tolerances)
- `update=true` — overwrite golden reference with current output (requires confirmation)
"""
function run_regression(case_dir; mode=:strict, update=false) end

"""
    run_all(; mode=:strict, filter=nothing)

Run all active cases. Optionally filter by name pattern.
Returns a RegressionReport.
"""
function run_all(; mode=:strict, filter=nothing) end

"""
    generate_reference(case_dir; generator=:julia)

Produce a new golden reference for a case.
Writes output to golden/ with a new version number and manifest.
"""
function generate_reference(case_dir; generator=:julia) end

"""
    compare(case_dir, version_a, version_b)

Compare two golden versions of the same case.
Useful for characterizing the impact of a physics change.
"""
function compare(case_dir, version_a, version_b) end

end # module
```


### Comparison Logic

The core comparison loads the golden NetCDF and the fresh run output, then checks each variable against its tolerance spec:

```julia
function compare_netcdf(golden_path, test_path, tolerances; verbose=true)
    results = Dict{String, ComparisonResult}()

    NCDataset(golden_path) do ds_gold
        NCDataset(test_path) do ds_test
            for varname in keys(tolerances)
                haskey(ds_gold, varname) || continue
                gold = ds_gold[varname][:]
                test = ds_test[varname][:]

                # Mask NaNs (both must be NaN in same locations)
                nan_match = isnan.(gold) .== isnan.(test)
                valid = .!isnan.(gold) .& .!isnan.(test)

                rtol = tolerances[varname].rtol
                atol = tolerances[varname].atol

                abs_err = abs.(gold[valid] .- test[valid])
                rel_err = abs_err ./ max.(abs.(gold[valid]), atol)

                max_abs = maximum(abs_err; init=0.0)
                max_rel = maximum(rel_err; init=0.0)
                mean_rel = length(rel_err) > 0 ? mean(rel_err) : 0.0
                pass = all(nan_match) && all(rel_err .<= rtol) && all(abs_err .<= atol)

                results[varname] = ComparisonResult(
                    pass = pass,
                    max_abs_error = max_abs,
                    max_rel_error = max_rel,
                    mean_rel_error = mean_rel,
                    nan_mismatch_count = count(.!nan_match),
                    n_compared = count(valid),
                    rtol = rtol,
                    atol = atol,
                )

                if verbose && !pass
                    @warn "FAIL: $varname" max_abs max_rel rtol atol
                    # Report worst-offending grid cells
                    worst_idx = argmax(rel_err)
                    # ... lat/lon of worst cell for debugging
                end
            end
        end
    end
    return results
end
```


### Streamline Comparison

Streamline CSV output is trickier because row counts can differ if the march terminates at slightly different steps. The comparison handles this:

```julia
function compare_streamlines(golden_csv, test_csv, tolerances)
    df_gold = CSV.read(golden_csv, DataFrame)
    df_test = CSV.read(test_csv, DataFrame)

    # Group by sink point (lat_sink, lon_sink, idx)
    # For each group:
    #   1. Check that both have the same number of steps (±1 allowed at truncation boundary)
    #   2. Compare tau, E0, dist, mu, wp at matched steps
    #   3. Allow the final 1-2 steps to differ (boundary effects from slightly different convergence)
end
```


### Golden Reference Versioning

When the model physics changes intentionally (e.g., you modify the τ̄ formula, change the interpolation scheme, or reactivate the EBM), the workflow is:

1. Run the regression suite in `science` mode to confirm changes are within expected bounds
2. If satisfied, bump the golden reference:
   ```
   julia> RegressionRunner.generate_reference("test/regression/cases/houston_bbox")
   # → writes golden/v002_CLIM.nc, golden/v002_STREAMLINES.csv, golden/v002_manifest.toml
   # → updates case.toml: current_version = "v002"
   ```
3. The manifest records the git commit, timestamp, generator (Python or Julia), and SHA-256 of the output files
4. Use `compare(case_dir, "v001", "v002")` to produce a diff report documenting what changed and by how much

Old golden versions are retained (but not checked against by default) so you can always trace the evolution of model output.


### CI Integration

```yaml
# .github/workflows/test.yml
name: Tests
on: [push, pull_request]
jobs:
  unit-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: julia-actions/setup-julia@v2
      - run: julia --project -e 'using Pkg; Pkg.test()'

  regression-strict:
    runs-on: ubuntu-latest
    steps:
      - uses: julia-actions/setup-julia@v2
      - run: |
          julia --project -e '
            include("test/regression/RegressionRunner.jl")
            report = RegressionRunner.run_all(mode=:strict)
            RegressionRunner.write_report(report, "test/regression/reports/ci_report.md")
            all(r.pass for r in values(report)) || exit(1)
          '
      - uses: actions/upload-artifact@v4
        if: failure()
        with:
          name: regression-report
          path: test/regression/reports/

  # Optional: GPU agreement test on self-hosted runner
  gpu-validation:
    runs-on: [self-hosted, gpu]
    if: github.ref == 'refs/heads/main'
    steps:
      - uses: julia-actions/setup-julia@v2
      - run: julia --project test/gpu/test_cpu_gpu_agreement.jl
```


### Adding a New Test Case (Checklist)

When you develop a new model configuration or expand to a new domain:

```
1. mkdir test/regression/cases/my_new_case/
2. Create case.toml (copy from houston_bbox, adjust parameters and tolerances)
3. Create generate_reference.jl (script that runs the model and saves output)
4. Run: julia generate_reference.jl → produces golden/v001_*.nc + manifest
5. Commit case.toml + generate_reference.jl (golden files via git-LFS or .gitignore)
6. The harness auto-discovers the new case on next run
```

No code changes to the test runner are needed — it scans `test/regression/cases/*/case.toml` at runtime.


### Dual-Tolerance Philosophy

The two tolerance modes serve different purposes through the project lifecycle:

- **`strict` mode** is for the port itself: Julia must reproduce Python's output to floating-point precision. This is the default in CI during Phase 1–2. Once the Julia code is validated, `strict` mode guards against accidental regressions from refactoring or optimization.

- **`science` mode** is for ongoing model development: when you intentionally change the physics (new interpolation scheme, modified τ̄ formula, EBM reactivation), strict tolerances will obviously break. Science mode lets you confirm the changes are within expected bounds before generating a new golden reference. The `case.toml` tolerance values document what "scientifically acceptable" means for each variable — a form of executable documentation.

---

## Estimated Performance Targets

| Configuration | Python (current) | Julia CPU (Phase 2) | Julia GPU (Phase 3) |
|--------------|-----------------|---------------------|---------------------|
| Small bbox (16×16, no decomp) | ~minutes | <1 second | <0.1 second |
| Regional (100×100, no decomp) | ~hours | ~1 minute | ~5 seconds |
| Regional (100×100, full decomp) | ~day | ~5 minutes | ~30 seconds |
| Global land (0.25°, no decomp) | infeasible | ~30 minutes | ~2 minutes |
| Global land (0.25°, full decomp) | infeasible | ~2 hours | ~10 minutes |

*Estimates assume M2 Ultra (24 CPU cores, 76 GPU cores) or NVIDIA A100. Actual performance depends on streamline convergence rates and memory bandwidth.*
