"""
    IO.jl — Data input/output.

Ported from Run.inputfiles() in attenuationMod_fxns.py, var_check() and
tau_arrayInitialize() from accessory_fxns.py, plus scattered xr I/O calls.
"""

using NCDatasets, DataFrames, CSV

# --- Data containers ---

"""
Container for climatology fields on a regular lat/lon grid.
All 2-D arrays are [lat, lon] (column-major friendly for meridional access).
"""
struct ClimatologyData{T<:AbstractFloat}
    lat::Vector{T}
    lon::Vector{T}
    P::Matrix{T}       # precipitation
    E::Matrix{T}       # evaporation
    UQ::Matrix{T}      # zonal moisture flux
    VQ::Matrix{T}      # meridional moisture flux
    Fmag::Matrix{T}    # moisture transport magnitude = sqrt(UQ^2 + VQ^2)
    LANDFRAC::Matrix{T}
end

"""
Container for 3-D climatology fields with a decomposition dimension.
Arrays are [lat, lon, decomp_dim].
"""
struct ClimatologyData3D{T<:AbstractFloat}
    lat::Vector{T}
    lon::Vector{T}
    P::Array{T,3}
    E::Array{T,3}
    UQ::Array{T,3}
    VQ::Array{T,3}
    Fmag::Array{T,3}
    LANDFRAC::Array{T,3}
end

# --- Variable checking ---

"""
    var_check(vars, ds)

Verify that all variable names in `vars` exist in the NCDataset `ds`.
Throws `ErrorException` if any are missing.
"""
function var_check(vars::Vector{String}, ds)
    for v in vars
        if !haskey(ds, v)
            error("Missing or mis-named variables. Required vars: $(join(vars, ", "))")
        end
    end
end

# --- Array initialization ---

"""
    tau_array_initialize(ds, lat, lon, config; decomp_dim=nothing)

Extract climatology arrays from NCDataset `ds`, matching
`tau_arrayInitialize` in accessory_fxns.py.

Returns `(P, E, UQ, VQ, Fmag, lfrac)` as plain arrays.
When `decomp_dim` is provided, arrays are 3-D [lat, lon, decomp].
"""
function tau_array_initialize(ds, lat::AbstractVector, lon::AbstractVector,
                              config::RunConfig;
                              decomp_dim::Union{Nothing,String}=nothing)
    e_name = config.vars.evaporation_field
    p_name = config.vars.precipitation_field
    u_name = config.vars.zonalqflux_field
    v_name = config.vars.meridqflux_field
    fmag_name = "Fmag"

    nlat = length(lat)
    nlon = length(lon)

    if decomp_dim === nothing
        P     = Array{Float64}(ds[p_name][:,:])
        E     = Array{Float64}(ds[e_name][:,:])
        UQ    = Array{Float64}(ds[u_name][:,:])
        VQ    = Array{Float64}(ds[v_name][:,:])
        Fmag  = Array{Float64}(ds[fmag_name][:,:])
        lfrac = Array{Float64}(ds["LANDFRAC"][:,:])
        return P, E, UQ, VQ, Fmag, lfrac
    else
        # 3-D: ensure [lat, lon, decomp] order
        P     = Array{Float64}(ds[p_name][:,:,:])
        E     = Array{Float64}(ds[e_name][:,:,:])
        UQ    = Array{Float64}(ds[u_name][:,:,:])
        VQ    = Array{Float64}(ds[v_name][:,:,:])
        Fmag  = Array{Float64}(ds[fmag_name][:,:,:])
        lfrac = Array{Float64}(ds["LANDFRAC"][:,:,:])
        return P, E, UQ, VQ, Fmag, lfrac
    end
end

# --- File I/O ---

"""
    load_climatology(config::RunConfig) -> NCDataset

Open the climatology NetCDF file. Caller is responsible for closing.
"""
function load_climatology(config::RunConfig)
    in_path = joinpath(config.io.run_path, config.io.run_name, config.io.input_dir)
    clim_path = joinpath(in_path, config.io.clim_fn)
    isfile(clim_path) || error("Cannot find climate input file: $clim_path")
    return NCDataset(clim_path, "r")
end

"""
    load_forcing(config::RunConfig) -> Union{NCDataset, Nothing}

Open the forcing NetCDF file, or return `nothing` if none configured.
"""
function load_forcing(config::RunConfig)
    config.io.force_fn === nothing && return nothing
    in_path = joinpath(config.io.run_path, config.io.run_name, config.io.input_dir)
    force_path = joinpath(in_path, config.io.force_fn)
    if isfile(force_path)
        return NCDataset(force_path, "r")
    else
        @warn "Could not find forcing file: $force_path — running isotope module only (no EBM)"
        return nothing
    end
end

"""
    load_coordinates(config::RunConfig) -> Union{DataFrame, Nothing}

Load coordinate list CSV when compute_method == "coord_list" or
collect_streamline_data == "some".
"""
function load_coordinates(config::RunConfig)
    iso = config.isotope
    needs_coords = (iso.compute_method == "coord_list") ||
                   (iso.collect_streamline_data == "some")
    needs_coords || return nothing

    in_path = joinpath(config.io.run_path, config.io.run_name, config.io.input_dir)
    coord_path = joinpath(in_path, iso.coord_list_filename)
    isfile(coord_path) || error("Cannot find coordinate list: $coord_path")

    df = CSV.read(coord_path, DataFrame)
    # Remove rows with missing lat/lon
    dropmissing!(df, [:lat, :lon])

    if iso.collect_streamline_data == "some"
        "collect_streamline" in names(df) ||
            error("Need 'collect_streamline' column in coord CSV when collecting 'some' streamlines")
    end

    return df
end

"""
    determine_runtype(config, ds_force) -> String

Determine run mode based on config and whether forcing data is available.
"""
function determine_runtype(config::RunConfig, ds_force)
    if ds_force !== nothing
        return config.isotope.solve_isotopes ? "full_model" : "EBM_noIsotopes"
    else
        return "isotopes_noEBM"
    end
end

"""
    save_results(outdata::Dict, streamlines::DataFrame, config::RunConfig)

Write output NetCDF and streamline CSV.
"""
function save_results(outdata::Dict{String,Any}, streamlines::DataFrame, config::RunConfig)
    results_dir = joinpath(config.io.run_path, config.io.run_name, "results")
    mkpath(results_dir)

    ds_path = joinpath(results_dir, config.io.run_name * "_CLIM.nc")
    strmln_path = joinpath(results_dir, config.io.run_name * "_STREAMLINES.csv")

    # Write NetCDF
    NCDataset(ds_path, "c") do ds
        lat = outdata["lat"]
        lon = outdata["lon"]
        defDim(ds, "lat", length(lat))
        defDim(ds, "lon", length(lon))
        ds_lat = defVar(ds, "lat", Float64, ("lat",))
        ds_lon = defVar(ds, "lon", Float64, ("lon",))
        ds_lat[:] = lat
        ds_lon[:] = lon

        for (name, arr) in outdata
            name in ("lat", "lon") && continue
            if ndims(arr) == 2
                v = defVar(ds, name, Float64, ("lat", "lon"))
                v[:,:] = arr
            elseif ndims(arr) == 3
                dim3_name = get(outdata, "dim3_name", "decomp")
                if !haskey(ds.dim, dim3_name)
                    defDim(ds, dim3_name, size(arr, 3))
                end
                v = defVar(ds, name, Float64, ("lat", "lon", dim3_name))
                v[:,:,:] = arr
            end
        end
    end

    # Write streamlines CSV
    CSV.write(strmln_path, streamlines)

    return ds_path, strmln_path
end
