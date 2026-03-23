#!/usr/bin/env julia
"""
    validate_vs_python.jl — Compare Julia iEBM output against Python reference.

Usage:
    julia --project scripts/validate_vs_python.jl <julia_nc> <python_nc> [--mode strict|science]
    julia --project scripts/validate_vs_python.jl --run-all <cases_dir>

Modes:
    strict  — port validation, rtol=1e-10 (default)
    science — physics evolution, rtol=1e-2
"""

# Add the regression module to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "test", "regression"))

using Printf
using Dates

include(joinpath(@__DIR__, "..", "test", "regression", "RegressionRunner.jl"))
using .RegressionRunner

function main()
    args = ARGS

    if isempty(args)
        println(stderr, "Usage: validate_vs_python.jl <julia_nc> <python_nc> [--mode strict|science]")
        println(stderr, "       validate_vs_python.jl --run-all <cases_dir>")
        exit(1)
    end

    if args[1] == "--run-all"
        cases_dir = length(args) >= 2 ? args[2] :
            joinpath(@__DIR__, "..", "test", "regression", "cases")

        println("Running all regression cases in: $cases_dir")
        results = run_all(cases_dir)

        if isempty(results)
            println("No cases found (or all skipped due to missing reference data).")
            println("Place reference data in the case directories to enable regression testing.")
            exit(0)
        end

        report_path = joinpath(@__DIR__, "..", "test", "regression", "reports",
                               "report_$(Dates.format(now(), "yyyymmdd_HHMMSS")).txt")
        report = generate_report(results; outpath=report_path)
        println(report)

        all_passed = all(r -> r.passed, results)
        exit(all_passed ? 0 : 1)
    end

    # Single file comparison mode
    if length(args) < 2
        println(stderr, "Error: need both julia and python output paths")
        exit(1)
    end

    julia_path = args[1]
    python_path = args[2]

    mode = :strict
    for i in 1:length(args)-1
        if args[i] == "--mode"
            mode = Symbol(args[i+1])
        end
    end

    if !(mode in (:strict, :science))
        println(stderr, "Error: mode must be 'strict' or 'science'")
        exit(1)
    end

    # Detect file type and compare
    if endswith(julia_path, ".nc")
        println("Comparing NetCDF files (mode=$mode):")
        println("  Julia:  $julia_path")
        println("  Python: $python_path")
        results = compare_netcdf(julia_path, python_path; mode=mode)
    elseif endswith(julia_path, ".csv")
        println("Comparing CSV files (mode=$mode):")
        println("  Julia:  $julia_path")
        println("  Python: $python_path")
        results = compare_streamlines(julia_path, python_path; mode=mode)
    else
        println(stderr, "Error: unsupported file type. Use .nc or .csv")
        exit(1)
    end

    # Print results
    println()
    println("=" ^ 60)
    all_passed = true
    for f in results
        status = f.passed ? "PASS" : "FAIL"
        @printf("  [%s] %-25s  max_rel=%.2e  max_abs=%.2e\n",
                status, f.name, f.max_rel_err, f.max_abs_err)
        if !f.passed
            all_passed = false
        end
    end
    println("=" ^ 60)
    println(all_passed ? "All fields match." : "SOME FIELDS FAILED — see above.")

    exit(all_passed ? 0 : 1)
end

main()
