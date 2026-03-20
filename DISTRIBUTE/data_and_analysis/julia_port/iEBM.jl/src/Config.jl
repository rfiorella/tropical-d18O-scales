"""
    Config.jl — Run configuration structs and defaults.

Ported from initialize.py (Houston case).
"""

# --- Sub-configs ---

Base.@kwdef struct GridConfig
    deg_per_lat::Float64 = 0.25
    deg_per_lon::Float64 = 0.25
end

Base.@kwdef struct VariableNames
    albedo_var::String            = "albedo_clr_surf"
    insolation_var::String        = "SOLIN"
    precip_input::String          = "PRECT"
    spechum_input::String         = "Q"
    evap_input::String            = "ET"
    zonalqflux_input::String      = "UQ"
    meridqflux_input::String      = "VQ"
    precipitation_field::String   = "PRECT"
    spechum_field::String         = "Q"
    evaporation_field::String     = "ET"
    zonalqflux_field::String      = "UQ"
    meridqflux_field::String      = "VQ"
    potentialevap_field::String   = "E0"
    budyko_omega_field::String    = "omega_budyko"
end

Base.@kwdef struct TopoConfig
    orog_partition::Bool              = false
    topo_varname::String              = "topo_m"
    topo_mask_varname::String         = "orog_mask"
    slope_threshold_m_per_cell::Float64 = 400.0
    elev_threshold_m::Float64         = 300.0
    use_ocean_in_orog_zonal_mean::Bool = true
end

Base.@kwdef struct ForcingConfig
    forcename_insol::String          = "insolation"
    force_multiplier_insol::Float64  = 1.0
    T_inertia_land::Float64          = 0.0
    T_inertia_ocean::Float64         = 1.0
    force_bbox_insol::Union{Nothing, Dict{String,Any}} = nothing
    forcename_albedo::String         = "albedo"
    force_multiplier_albedo::Float64 = 1.0
    TOAalbedo_sensitivity_factor::Float64 = 0.7
    forcename_arbitrary::String      = "forcing"
    force_multiplier_arbitrary::Float64 = 1.0
    arbitrary_forcing_Wm2::Float64   = 0.0
    arbitrary_forcing_landOnly::Bool = true
    force_bbox_arbitrary::Union{Nothing, Dict{String,Any}} = nothing
end

Base.@kwdef struct EBMConfig
    EFE_absmax::Float64             = 35.0
    heaviside_deglim_V::Float64     = 20.0
    EFE_threshold::Float64          = 0.1e8
    EFE_threshold_step::Float64     = 0.02e8
    EFE_threshold_max::Float64      = 0.5e8
    DualEFE_ChiRangeForMean::Float64 = 0.2e8
    EFPM_nLinesMax_CTRL::Int        = 2
    EFPM_nLinesMax_MATCH::Int       = 3
    EFPM_lat_absmax::Float64        = 45.0
    EFPM_stitch_search_radius_deg::Float64 = 17.0
    EFPMs_Match_ProximityLimit::Float64 = 50.0
    EFPM_removeHorizontalLines::Bool = true
    EFPM_HorizontalLine_angle::Float64 = 30.0
    EFPM_ChiThresholdFactor::Float64 = 0.2
    EFPM_MinLatSpan::Float64        = 20.0
    heaviside_deglim_U::Float64     = 70.0
    lat_limit_zonal_P_shift::Float64 = 45.0
end

Base.@kwdef struct IsotopeConfig
    solve_isotopes::Bool                    = true
    compute_method::String                  = "bbox"   # "bbox" or "coord_list"
    coord_list_filename::String             = "scratch_streamlines2save.csv"
    compute_land_only::Bool                 = true
    bbox_lat_range::Vector{Float64}         = [28.0, 32.0]
    bbox_lon_range::Vector{Float64}         = [263.0, 267.0]
    bbox_resolution::Vector{Float64}        = [0.25, 0.25]
    streamline_dx_km::Float64               = 14.0
    streamline_max_tau::Float64             = 8.0
    streamline_max_dist_km::Float64         = 25000.0
    n_samples_taubar_land::Int              = 10
    collect_streamline_data::String         = "some"   # "all", "some", "none"
    streamline_save_coarsener::Int          = 10
end

Base.@kwdef struct DecompConfig
    local_v_regional_local_evap::Bool       = false
    local_v_regional_upwind::Bool           = false
    local_threshold_km::Float64             = 1000.0
    solve_frac_dist_local_evap::Bool        = false
    solve_frac_dist_streamline::Bool        = false
    decomp_E_L_s::Bool                      = false
    ELs_initstate_same_yrslice::Bool        = false
    dtau_fraction::Float64                  = 0.75
    save_local_v_regional_streamlines::Bool = false
end

Base.@kwdef struct IOConfig
    run_name::String       = "houston"
    run_path::String       = "."
    input_dir::String      = "input/"
    clim_fn::String        = "era_mon_fixvars.nc"
    force_fn::Union{Nothing,String} = nothing
    slices_to_solve::Union{Nothing,Vector{String}} = nothing
    filesave_suffix::String = ""
end

Base.@kwdef struct PhysicalConstants
    Lv::Float64 = 2.54e6    # [J kg-1] latent heat of vaporization
    g::Float64  = 9.81      # [m s-2]  gravitational acceleration
end

# --- Top-level config ---

Base.@kwdef struct RunConfig
    io::IOConfig               = IOConfig()
    grid::GridConfig           = GridConfig()
    vars::VariableNames        = VariableNames()
    topo::TopoConfig           = TopoConfig()
    forcing::ForcingConfig     = ForcingConfig()
    ebm::EBMConfig             = EBMConfig()
    isotope::IsotopeConfig     = IsotopeConfig()
    decomp::DecompConfig       = DecompConfig()
    constants::PhysicalConstants = PhysicalConstants()
end

"""
    validate(config::RunConfig)

Check parameter consistency; throw on invalid combinations.
"""
function validate(config::RunConfig)
    iso = config.isotope
    dc  = config.decomp

    iso.compute_method in ("bbox", "coord_list") ||
        error("compute_method must be \"bbox\" or \"coord_list\", got \"$(iso.compute_method)\"")
    iso.collect_streamline_data in ("all", "some", "none") ||
        error("collect_streamline_data must be \"all\", \"some\", or \"none\"")

    if length(iso.bbox_lat_range) != 2
        error("bbox_lat_range must have exactly 2 elements")
    end
    if length(iso.bbox_lon_range) != 2
        error("bbox_lon_range must have exactly 2 elements")
    end

    iso.streamline_dx_km > 0 || error("streamline_dx_km must be > 0")
    iso.streamline_max_tau > 0 || error("streamline_max_tau must be > 0")
    iso.streamline_max_dist_km > 0 || error("streamline_max_dist_km must be > 0")
    iso.n_samples_taubar_land >= 1 || error("n_samples_taubar_land must be >= 1")

    0.0 <= dc.dtau_fraction <= 1.0 || error("dtau_fraction must be in [0, 1]")

    return nothing
end

"""
    make_houston_config() -> RunConfig

Convenience: return the default Houston case configuration.
"""
make_houston_config() = RunConfig()
