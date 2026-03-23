#!/usr/bin/env julia
"""
    bench_fullrun.jl — Benchmark a full iEBM run.

Requires input data to be available. Skips gracefully if not found.
"""

using Printf

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using iEBM

function main()
    config = make_houston_config()

    # Check if input data exists
    input_file = config.io.input_file
    if !isfile(input_file)
        println("Input data not found at: $input_file")
        println("Skipping full-run benchmark. Provide data to enable.")
        return
    end

    println("=" ^ 50)
    println("  iEBM Full Run Benchmark")
    println("  Input: $input_file")
    println("=" ^ 50)

    # Warmup
    println("\nWarmup run...")
    t_warmup = @elapsed runloop(config)
    @printf("  Warmup: %.2f s\n", t_warmup)

    # Timed runs
    n_trials = 3
    times = Float64[]
    for i in 1:n_trials
        println("Trial $i / $n_trials ...")
        t = @elapsed runloop(config)
        push!(times, t)
        @printf("  %.2f s\n", t)
    end

    println()
    @printf("Results (n=%d):\n", n_trials)
    @printf("  min:  %.2f s\n", minimum(times))
    @printf("  mean: %.2f s\n", sum(times) / n_trials)
    @printf("  max:  %.2f s\n", maximum(times))
end

main()
