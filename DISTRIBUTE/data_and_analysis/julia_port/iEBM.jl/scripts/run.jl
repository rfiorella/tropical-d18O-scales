#!/usr/bin/env julia
"""
    run.jl — Entry point for running the iEBM model.

Replaces RUN.py. Configure via RunConfig and call runloop().

Usage:
    julia --project=iEBM.jl scripts/run.jl
"""

using iEBM

# Build configuration (Houston case defaults)
config = RunConfig(
    io = IOConfig(
        run_name   = "houston",
        run_path   = joinpath(@__DIR__, "..", "..", ".."),
        input_dir  = "input/",
        clim_fn    = "era_mon_fixvars.nc",
        force_fn   = nothing,
    ),
    isotope = IsotopeConfig(
        compute_method      = "bbox",
        bbox_lat_range      = [28.0, 32.0],
        bbox_lon_range      = [263.0, 267.0],
        bbox_resolution     = [0.25, 0.25],
        streamline_dx_km    = 14.0,
        streamline_max_tau  = 8.0,
        streamline_max_dist_km = 25000.0,
        n_samples_taubar_land  = 10,
        collect_streamline_data = "some",
        streamline_save_coarsener = 10,
    ),
)

# Run
output, streamlines = runloop(config)

println("Output fields: ", keys(output))
println("Streamline rows: ", nrow(streamlines))
