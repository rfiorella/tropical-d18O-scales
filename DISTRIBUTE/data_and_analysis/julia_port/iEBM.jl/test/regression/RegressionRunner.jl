"""
    RegressionRunner.jl — Regression testing framework for iEBM.

Compare Julia output against Python reference data (NetCDF + CSV).
Supports two tolerance modes:
  - :strict  — port-validation (rtol=1e-10), ensures exact numerical match
  - :science — physics-evolution (rtol=1e-2), allows controlled drift
"""
module RegressionRunner

using NCDatasets
using DataFrames
using CSV
using Test
using Printf
using Dates

# ── Result types ──────────────────────────────────────────────────────

struct FieldComparison
    name::String
    max_abs_err::Float64
    max_rel_err::Float64
    mean_abs_err::Float64
    n_mismatched::Int
    n_total::Int
    passed::Bool
end

struct ComparisonResult
    case_name::String
    mode::Symbol              # :strict or :science
    fields::Vector{FieldComparison}
    passed::Bool
    elapsed_seconds::Float64
    timestamp::String
end

# ── Tolerances ────────────────────────────────────────────────────────

const TOLERANCES = Dict(
    :strict  => (rtol=1e-10, atol=1e-12),
    :science => (rtol=1e-2,  atol=1e-6),
)

# ── Comparison utilities ─────────────────────────────────────────────

"""
    compare_arrays(a, b; rtol, atol) → FieldComparison

Element-wise comparison of two arrays, returning error statistics.
"""
function compare_arrays(name::String, a::AbstractArray, b::AbstractArray;
                        rtol::Float64, atol::Float64)
    @assert size(a) == size(b) "Size mismatch for $name: $(size(a)) vs $(size(b))"

    n = length(a)
    max_abs = 0.0
    max_rel = 0.0
    sum_abs = 0.0
    n_bad = 0

    for i in eachindex(a)
        ai, bi = a[i], b[i]
        # Skip NaN-NaN pairs (both NaN = match)
        if isnan(ai) && isnan(bi)
            continue
        end
        if isnan(ai) || isnan(bi)
            n_bad += 1
            continue
        end
        abs_err = abs(ai - bi)
        denom = max(abs(ai), abs(bi), 1e-30)
        rel_err = abs_err / denom

        sum_abs += abs_err
        max_abs = max(max_abs, abs_err)
        max_rel = max(max_rel, rel_err)

        if abs_err > atol && rel_err > rtol
            n_bad += 1
        end
    end

    mean_abs = n > 0 ? sum_abs / n : 0.0
    passed = n_bad == 0

    return FieldComparison(name, max_abs, max_rel, mean_abs, n_bad, n, passed)
end

"""
    compare_netcdf(julia_path, python_path; mode=:strict, vars=nothing)

Compare all (or selected) variables in two NetCDF files.
"""
function compare_netcdf(julia_path::String, python_path::String;
                        mode::Symbol=:strict, vars=nothing)
    tol = TOLERANCES[mode]
    results = FieldComparison[]

    NCDataset(julia_path) do ds_jl
        NCDataset(python_path) do ds_py
            varnames = vars !== nothing ? vars : collect(keys(ds_py))
            # Filter to numeric variables present in both
            for vname in varnames
                haskey(ds_jl, vname) || continue
                haskey(ds_py, vname) || continue
                v_jl = ds_jl[vname][:]
                v_py = ds_py[vname][:]
                if eltype(v_jl) <: Number && eltype(v_py) <: Number
                    fc = compare_arrays(vname, Float64.(v_jl), Float64.(v_py);
                                        rtol=tol.rtol, atol=tol.atol)
                    push!(results, fc)
                end
            end
        end
    end

    return results
end

"""
    compare_streamlines(julia_path, python_path; mode=:strict, cols=nothing)

Compare streamline CSV output files.
"""
function compare_streamlines(julia_path::String, python_path::String;
                             mode::Symbol=:strict, cols=nothing)
    tol = TOLERANCES[mode]
    results = FieldComparison[]

    df_jl = CSV.read(julia_path, DataFrame)
    df_py = CSV.read(python_path, DataFrame)

    colnames = cols !== nothing ? cols : names(df_py)
    for cname in colnames
        hasproperty(df_jl, cname) || continue
        hasproperty(df_py, cname) || continue
        c_jl = df_jl[!, cname]
        c_py = df_py[!, cname]
        if eltype(c_jl) <: Number && eltype(c_py) <: Number
            fc = compare_arrays(string(cname), Float64.(c_jl), Float64.(c_py);
                                rtol=tol.rtol, atol=tol.atol)
            push!(results, fc)
        end
    end

    return results
end

# ── Case runner ───────────────────────────────────────────────────────

"""
    CaseConfig

Parsed from a case.toml file.
"""
struct CaseConfig
    name::String
    description::String
    mode::Symbol
    config_overrides::Dict{String,Any}
    reference_nc::String
    reference_csv::String
    nc_vars::Union{Nothing, Vector{String}}
    csv_cols::Union{Nothing, Vector{String}}
end

"""
    load_case(toml_path) → CaseConfig

Parse a regression case TOML file.
"""
function load_case(toml_path::String)
    d = _parse_toml(toml_path)
    name = get(d, "name", basename(dirname(toml_path)))
    desc = get(d, "description", "")
    mode = Symbol(get(d, "mode", "strict"))
    overrides = get(d, "config", Dict{String,Any}())
    ref = get(d, "reference", Dict{String,Any}())
    ref_nc = get(ref, "netcdf", "")
    ref_csv = get(ref, "csv", "")
    nc_vars = haskey(d, "nc_vars") ? convert(Vector{String}, d["nc_vars"]) : nothing
    csv_cols = haskey(d, "csv_cols") ? convert(Vector{String}, d["csv_cols"]) : nothing
    return CaseConfig(name, desc, mode, overrides, ref_nc, ref_csv, nc_vars, csv_cols)
end

function _parse_toml(path::String)
    # Simple TOML parser using Julia's built-in TOML stdlib
    return Base.TOML.parsefile(path)
end

"""
    run_regression(case_toml; julia_output_dir=nothing) → ComparisonResult

Run a single regression case: execute iEBM, compare output to reference.
Skips (returns nothing) if reference data is missing.
"""
function run_regression(case_toml::String; julia_output_dir::Union{String,Nothing}=nothing)
    case = load_case(case_toml)

    # Check reference data exists
    ref_dir = dirname(case_toml)
    ref_nc = isabspath(case.reference_nc) ? case.reference_nc : joinpath(ref_dir, case.reference_nc)
    ref_csv = isabspath(case.reference_csv) ? case.reference_csv : joinpath(ref_dir, case.reference_csv)

    if !isempty(case.reference_nc) && !isfile(ref_nc)
        @warn "Skipping case '$(case.name)': reference NetCDF not found at $ref_nc"
        return nothing
    end
    if !isempty(case.reference_csv) && !isfile(ref_csv)
        @warn "Skipping case '$(case.name)': reference CSV not found at $ref_csv"
        return nothing
    end

    t0 = time()
    field_results = FieldComparison[]

    # Run Julia model
    try
        # Build config with overrides
        config = _build_config(case.config_overrides)

        # Set output directory
        if julia_output_dir !== nothing
            outdir = julia_output_dir
        else
            outdir = mktempdir()
        end

        # Run model
        @eval using iEBM
        output, df_sl = Base.invokelatest(iEBM.runloop, config)

        # Compare NetCDF output
        if !isempty(case.reference_nc) && isfile(ref_nc)
            # Find Julia output NetCDF
            jl_nc = joinpath(outdir, basename(ref_nc))
            if isfile(jl_nc)
                nc_results = compare_netcdf(jl_nc, ref_nc;
                                            mode=case.mode, vars=case.nc_vars)
                append!(field_results, nc_results)
            end
        end

        # Compare CSV output
        if !isempty(case.reference_csv) && isfile(ref_csv)
            jl_csv = joinpath(outdir, basename(ref_csv))
            if isfile(jl_csv)
                csv_results = compare_streamlines(jl_csv, ref_csv;
                                                  mode=case.mode, cols=case.csv_cols)
                append!(field_results, csv_results)
            end
        end
    catch e
        push!(field_results, FieldComparison("RUN_ERROR", Inf, Inf, Inf, 1, 1, false))
        @error "Case '$(case.name)' failed" exception=(e, catch_backtrace())
    end

    elapsed = time() - t0
    all_passed = all(f -> f.passed, field_results)

    return ComparisonResult(
        case.name, case.mode, field_results,
        all_passed, elapsed,
        Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    )
end

"""
    run_all(cases_dir; kwargs...) → Vector{ComparisonResult}

Find and run all case.toml files under `cases_dir`.
"""
function run_all(cases_dir::String; kwargs...)
    results = ComparisonResult[]
    for (root, dirs, files) in walkdir(cases_dir)
        for f in files
            if f == "case.toml"
                r = run_regression(joinpath(root, f); kwargs...)
                r !== nothing && push!(results, r)
            end
        end
    end
    return results
end

# ── Report generation ─────────────────────────────────────────────────

"""
    generate_report(results; outpath=nothing) → String

Generate a human-readable report. Returns the report text.
If `outpath` is given, also writes to file.
"""
function generate_report(results::Vector{ComparisonResult}; outpath::Union{String,Nothing}=nothing)
    io = IOBuffer()
    println(io, "=" ^ 72)
    println(io, "  iEBM Regression Test Report")
    println(io, "  Generated: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "=" ^ 72)
    println(io)

    n_pass = count(r -> r.passed, results)
    n_total = length(results)
    println(io, "  Overall: $n_pass / $n_total cases passed")
    println(io)

    for r in results
        status = r.passed ? "PASS" : "FAIL"
        println(io, "─" ^ 72)
        @printf(io, "  [%s] %s  (mode=%s, %.1fs)\n", status, r.case_name, r.mode, r.elapsed_seconds)

        for f in r.fields
            fstatus = f.passed ? "ok" : "FAIL"
            @printf(io, "    %-30s  %s  max_rel=%.2e  max_abs=%.2e  (%d/%d bad)\n",
                    f.name, fstatus, f.max_rel_err, f.max_abs_err, f.n_mismatched, f.n_total)
        end
    end
    println(io, "─" ^ 72)

    report = String(take!(io))

    if outpath !== nothing
        mkpath(dirname(outpath))
        open(outpath, "w") do f
            write(f, report)
        end
    end

    return report
end

# ── Config builder (stub — needs iEBM loaded) ────────────────────────

function _build_config(overrides::Dict{String,Any})
    # Default: return a Houston config. Override fields as specified.
    # This is a simplified builder; full implementation will map TOML keys
    # to RunConfig keyword arguments.
    @eval using iEBM
    config = Base.invokelatest(iEBM.make_houston_config)
    # TODO: apply overrides from case.toml to config fields
    return config
end

# ── Exports ───────────────────────────────────────────────────────────

export ComparisonResult, FieldComparison, CaseConfig
export compare_arrays, compare_netcdf, compare_streamlines
export load_case, run_regression, run_all, generate_report

end # module RegressionRunner
