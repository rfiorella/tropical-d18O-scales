#!/usr/bin/env julia
"""
    bench_streamline.jl — Benchmark the streamline marching kernel.

Measures march_streamline! performance on synthetic uniform flow fields
at various grid sizes. No input data required.
"""

using Printf

# Activate the project
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using iEBM

function bench_march(; Nmax=500, n_warmup=5, n_trials=100)
    # Build a uniform flow field for benchmarking
    lat = collect(Float64, -30:1.0:30)
    lon = collect(Float64, 0:1.0:359)
    nlat = length(lat)
    nlon = length(lon)

    P = fill(1e-5, nlat, nlon)
    E = fill(1e-5, nlat, nlon)
    UQ = fill(10.0, nlat, nlon)
    VQ = fill(0.001, nlat, nlon)
    lfrac = ones(nlat, nlon)

    Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = build_interpolators(
        lat, lon, E, P, UQ, VQ, lfrac)

    dx = 14.0
    Dx = dx / 111.0
    result = StreamlineResult(Nmax)

    # Warmup
    for _ in 1:n_warmup
        march_streamline!(result, 0.0, 180.0, Nmax, 8.0, Dx, dx,
                          Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)
    end

    # Timed runs
    times = Float64[]
    for _ in 1:n_trials
        t = @elapsed march_streamline!(result, 0.0, 180.0, Nmax, 8.0, Dx, dx,
                                        Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)
        push!(times, t)
    end

    med = sort(times)[div(n_trials, 2)]
    mn = minimum(times)
    avg = sum(times) / n_trials

    @printf("march_streamline! (Nmax=%d)\n", Nmax)
    @printf("  min:    %.2f μs\n", mn * 1e6)
    @printf("  median: %.2f μs\n", med * 1e6)
    @printf("  mean:   %.2f μs\n", avg * 1e6)
    @printf("  steps:  %d\n", result.nsteps)

    return (min=mn, median=med, mean=avg)
end

function bench_interpolator(; n_warmup=100, n_trials=10_000)
    lat = collect(Float64, -90:1.0:90)
    lon = collect(Float64, 0:1.0:359)
    nlat = length(lat)
    nlon = length(lon)

    data = [sin(deg2rad(la)) * cos(deg2rad(lo)) for la in lat, lo in lon]
    P = fill(1e-5, nlat, nlon)
    UQ = fill(10.0, nlat, nlon)
    VQ = fill(0.001, nlat, nlon)
    lfrac = ones(nlat, nlon)

    Efit, _, _, _, _, _ = build_interpolators(lat, lon, data, P, UQ, VQ, lfrac)

    # Warmup
    for _ in 1:n_warmup
        Efit(15.3, 42.7)
    end

    # Timed runs
    times = Float64[]
    for _ in 1:n_trials
        t = @elapsed Efit(15.3, 42.7)
        push!(times, t)
    end

    med = sort(times)[div(n_trials, 2)]
    mn = minimum(times)

    @printf("RegularGridInterp call\n")
    @printf("  min:    %.0f ns\n", mn * 1e9)
    @printf("  median: %.0f ns\n", med * 1e9)

    return (min=mn, median=med)
end

function main()
    println("=" ^ 50)
    println("  iEBM Streamline Benchmark")
    println("=" ^ 50)
    println()

    bench_interpolator()
    println()

    for nmax in [100, 200, 500, 1000]
        bench_march(Nmax=nmax)
        println()
    end
end

main()
