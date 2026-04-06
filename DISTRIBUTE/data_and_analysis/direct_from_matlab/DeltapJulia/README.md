# DeltapModel.jl

Julia port of `get_deltap_fast_optimized.m` -- a Lagrangian streamline-based
Rayleigh distillation model for precipitation delta-18O. Computes `deltap_bar`
(precipitation isotope ratio) and `tau_bar` (residence-time metric) on a
regular lat/lon grid.

Three backends are provided:

| Function | Description |
|---|---|
| `get_deltap_cpu` | Naive single-threaded CPU (reference implementation) |
| `get_deltap_cpu_threaded` | O(1) interpolation + `@threads` (recommended) |
| `run_gpu` | KernelAbstractions.jl kernel (Metal / CUDA / CPU fallback) |

## Requirements

- **Julia >= 1.12**
- **GNU Octave** (only needed to regenerate reference fixtures from scratch)

All Julia dependencies are declared in `Project.toml` and installed
automatically by the package manager.

## Quick start

```bash
cd DeltapJulia/

# Install dependencies (first time only)
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run all tests
julia --project -e 'using Pkg; Pkg.test()'
```

## Project layout

```
DeltapJulia/
├── Project.toml                  # Package manifest
├── Manifest.toml                 # Locked dependency versions
├── data/
│   └── alpha_eq.csv              # Equilibrium fractionation factors
├── src/
│   ├── DeltapModel.jl            # Module entry point
│   ├── interpolation.jl          # Bilinear interp + RegularGrid (O(1))
│   ├── inpaint.jl                # NaN infill (iterative Laplacian)
│   ├── types.jl                  # DeltapFields struct for GPU
│   ├── streamline_cpu.jl         # get_deltap_cpu (naive)
│   ├── streamline_threaded.jl    # get_deltap_cpu_threaded
│   ├── streamline_gpu.jl         # @kernel deltap_kernel! (two-trace)
│   └── gpu_driver.jl             # run_gpu + backend selection
├── test/
│   ├── runtests.jl               # Test orchestrator
│   ├── test_interpolation.jl     # Unit tests for interpolation
│   ├── test_reference.jl         # Naive CPU vs Octave reference
│   ├── test_reference_threaded.jl# Threaded CPU vs Octave reference
│   ├── test_reference_gpu.jl     # GPU kernel vs Octave reference
│   └── fixtures/                 # Reference .mat files (see below)
├── bench/
│   ├── bench_cpu.jl              # Small-grid benchmark
│   └── bench_halfdeg.jl          # 0.5-degree profiling script
└── scripts/
    ├── generate_synthetic_testcase.m   # Octave: tiny + small inputs
    ├── generate_reference.m            # Octave: run model, save outputs
    ├── generate_halfdeg_testcase.m     # Octave: 0.5-degree inputs
    ├── generate_halfdeg_reference.jl   # Julia: 0.5-degree reference
    ├── get_deltap_fast_optimized_octave.m  # Octave-compat model
    └── inpaint_nans.m                  # Octave-compat NaN infill
```

## Tests

### Running the full test suite

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

This runs **88 tests** covering:

| Test file | What it checks | Count |
|---|---|---|
| `test_interpolation.jl` | Bilinear interp (searchsorted + RegularGrid), 1D interp, RegularLookup, edge clamping | 22 |
| `test_reference.jl` | Naive CPU vs Octave `.mat` reference (tiny + small + halfdeg) | 18 |
| `test_reference_threaded.jl` | Threaded CPU vs Octave reference | 18 |
| `test_reference_gpu.jl` | GPU kernel (CPU backend) vs reference, both Float64 and Float32 | 30 |

**Tolerance thresholds:**

- Float64 vs Octave reference: max absolute error < 1e-4 (typical: ~5e-7)
- Float32 (GPU) vs reference: max absolute error < 1e-3 (typical: ~5e-4)

### Reference fixtures

The test suite requires `.mat` files in `test/fixtures/`. These are
pre-generated and checked into the repo:

| File | Grid | Points | Source |
|---|---|---|---|
| `reference_tiny.mat` | 12x6 | 72 | Octave |
| `reference_small.mat` | 60x30 | 1,800 | Octave |
| `reference_halfdeg.mat` | 720x360 | 259,200 | Julia naive CPU |

Each file contains all model inputs (`E`, `P`, `UQ`, `VQ`, `Tcond`, `LAT`,
`LON`, `LAT2`, `LON2`, `delta_e`, `alpha_eq`, `Plim`, `dmax`) plus the
reference outputs (`deltap_bar`, `tau_bar`).

### Regenerating fixtures from scratch

If you need to regenerate the reference data:

**Tiny and small cases (requires GNU Octave):**

```bash
cd scripts/

# Generate synthetic input fields
octave --no-gui generate_synthetic_testcase.m

# Run the Octave model to produce reference outputs
octave --no-gui generate_reference.m
```

**Half-degree case (Julia, much faster than Octave):**

```bash
cd scripts/

# Generate synthetic inputs (Octave)
octave --no-gui generate_halfdeg_testcase.m

# Generate reference outputs (Julia naive CPU, ~12s)
julia --project=.. generate_halfdeg_reference.jl
```

### Running a single test file

```bash
# Just the interpolation unit tests
julia --project -e 'using Test; include("test/test_interpolation.jl")'

# Just the GPU reference tests
julia --project -e '
    using Test, MAT, DeltapModel, KernelAbstractions
    include("test/test_reference_gpu.jl")
'
```

## Benchmarks

### Quick benchmark (tiny + small grids)

```bash
# Single-threaded
julia --project bench/bench_cpu.jl

# With 4 threads
julia --project -t4 bench/bench_cpu.jl
```

### Full profiling (0.5-degree global grid)

```bash
# Single-threaded
julia --project bench/bench_halfdeg.jl

# Multi-threaded (recommended)
julia --project -t8 bench/bench_halfdeg.jl
```

The halfdeg benchmark tests all backends (naive CPU, threaded CPU, GPU kernel
on CPU, Metal GPU) and tunes workgroup sizes. It also validates each backend
against the reference data.

### Typical results (Apple M3, 259,200 grid points)

| Backend | Time | Speedup |
|---|---|---|
| Naive CPU (1 thread) | 11.2s | 1x |
| Threaded CPU (1 thread) | 4.4s | 2.6x |
| Threaded CPU (4 threads) | 1.6s | 7x |
| Threaded CPU (8 threads) | 1.1s | 10x |
| Metal GPU (Float32) | 73s | 0.15x |

The **threaded CPU backend is recommended** for production use. The GPU kernel
is slower for this workload due to variable-length streamlines causing SIMD
divergence and scattered memory access patterns from trajectory-following
interpolation.

## Usage

```julia
using DeltapModel

# Load your data (all 2D fields must have size (nlon, nlat))
# ...

# Recommended: threaded CPU
deltap_bar, tau_bar = get_deltap_cpu_threaded(
    E, P, UQ, VQ, Tcond,
    LAT, LON, LAT2, LON2,
    delta_e, alpha_eq, Plim, dmax
)

# GPU (auto-selects Metal if available, else CPU fallback)
deltap_bar, tau_bar = run_gpu(
    E, P, UQ, VQ, Tcond,
    LAT, LON, LAT2, LON2,
    delta_e, alpha_eq, Plim, dmax;
    T=Float32
)
```

**Input convention (matching Matlab):** all 2D fields have size `(nlon, nlat)`
where dimension 1 (rows) = longitude, dimension 2 (columns) = latitude.

**`alpha_eq`** is an Nx2 matrix: column 1 = temperature (degrees C), column 2 =
equilibrium fractionation factor. A CSV file is provided in `data/alpha_eq.csv`.

**`Plim`** is a 2D precipitation threshold array (same size as `P`). Grid
points where `P < Plim` are skipped (output is NaN).

**`dmax`** is the maximum upstream integration distance in km (e.g., 3000).

## Grid convention note

The Matlab source uses `[LAT, LON] = meshgrid(lat1d, lon1d)` which produces
arrays where rows index longitude and columns index latitude. This convention
is preserved in the Julia port. Internally, fields are transposed to
`(nlat, nlon)` for interpolation, with cyclic extension along the longitude
dimension.
