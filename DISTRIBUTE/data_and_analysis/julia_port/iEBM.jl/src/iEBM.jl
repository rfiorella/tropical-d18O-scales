"""
    iEBM.jl — Isotope-enabled Energy Balance / Attenuation Model

Main module. Includes all submodules and defines public API.
"""
module iEBM

# Core dependencies
using NCDatasets
using DataFrames
using CSV

# Module files (order matters: dependencies first)
include("Config.jl")
include("Interpolation.jl")
include("Streamlines.jl")
include("IO.jl")
include("Hydroclim.jl")
include("EBM.jl")
include("TauBar.jl")
include("Decomposition.jl")
include("GPUKernels.jl")

# --- Exports ---

# Config
export RunConfig, GridConfig, VariableNames, TopoConfig, ForcingConfig,
       EBMConfig, IsotopeConfig, DecompConfig, IOConfig, PhysicalConstants
export validate, make_houston_config

# IO
export ClimatologyData, ClimatologyData3D
export var_check, tau_array_initialize
export load_climatology, load_forcing, load_coordinates, determine_runtype, save_results

# Interpolation
export RegularGridInterp, build_interpolators

# Hydroclim
export compute_slope, build_orog_mask, partition_field

# EBM stubs
export solve_ebm!, compute_efe!, compute_efpm!

# Streamlines
export StreamlineResult, march_streamline!, reset!

# TauBar
export linear_idx_to_degrees, dist_to_tau_coords, terrestrial_efrac
export compute_taubar_at_point, compute_taubar!, build_streamlines_df

# Decomposition
export spatial_decompose, local_plus_regional_streamline
export find_fraction_distance_localevap, find_fraction_distance_streamline
export ELpath_decompose, decompose_taubar!

# GPU stubs
export CPUBackend, select_backend, to_device

# --- Run orchestration ---

"""
    runloop(config::RunConfig)

Main entry point: load data, run model, save results.
Port of Run.runloop() from attenuationMod_fxns.py.
"""
function runloop(config::RunConfig)
    validate(config)

    println("********************************************")
    println("Now starting run: $(config.io.run_name)")
    println("--------------------------------------------")

    # Load inputs
    ds_clim = load_climatology(config)
    ds_force = load_forcing(config)
    df_coords = load_coordinates(config)
    runtype = determine_runtype(config, ds_force)

    dc = config.decomp
    spatial_decomposition = dc.local_v_regional_local_evap || dc.local_v_regional_upwind
    clim_decomposition = dc.decomp_E_L_s
    tau_decompose = spatial_decomposition || clim_decomposition

    if !tau_decompose
        output, df_streamlines = compute_taubar!(ds_clim, config;
                                                  iso_coords=df_coords,
                                                  runtype=runtype)
    else
        output, df_streamlines = decompose_taubar!(ds_clim, config;
                                                    iso_coords=df_coords,
                                                    runtype=runtype)
    end

    # Close datasets
    close(ds_clim)
    ds_force !== nothing && close(ds_force)

    # Save results
    save_results(output, df_streamlines, config)

    println("********************************************")
    println("run complete!                           :)")
    println("********************************************")

    return output, df_streamlines
end

export runloop

end # module
