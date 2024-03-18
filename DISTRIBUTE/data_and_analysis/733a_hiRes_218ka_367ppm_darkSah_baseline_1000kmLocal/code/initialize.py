# %%
# ------------------------------------------- # 
# RUN 
# ------------------------------------------- #

# %% 

init = {

    # CASE SETP + INPUTS --------------------------------- #
    # ---------------------------------------------------- #
    "run_name":   "733a_hiRes_218ka_367ppm_darkSah_baseline_1000kmLocal",                # must also be dir name
    "run_path":   "/Users/tylerkukla/Documents/academic_research/myProjects/spatialScales_tropicalIsotopes/DISTRIBUTE/data_and_analysis", # where run_name lies
    "input_dir":  "inputs/",                                 # where input files lie (w/in run_name)
    "clim_fn":    "mio_darkSah_218ka_367ppm+1850alb_367ppm_mjjas_ndjfm_input.nc",          # input climatology
    "force_fn":   None,                  # input forcing
    "slices_to_solve":  None,             # names of yr_slices to solve (will solve all slices if empty)               

    # FILE-FORMATTING ------------------------------------ # 
    # ---------------------------------------------------- #
    'filesave_suffix': "",



    # CLIMATOLOGY ---------------------------------------- #
    # ---------------------------------------------------- #

    # --- GRID --------------
    "deg_per_lat": 1.875, 
    "deg_per_lon": 2.5,

    # --- VARIABLES ---------
    "albedo_var": "albedo_clr_surf",   # [] name of albedo variable
    "insolation_var": "SOLIN",         # [] name of insolation variable
    "Inputfile_name_Precip": "PRECT",  # [] name of precipitation var
    "Inputfile_name_SpecHum": "Q",     # [] name of specific humidity var
    "Inputfile_name_Evap": "ET",       # [] name of evaporation var
    # ... these names may update throughout model run
    "Precipitation_Field": "PRECT",      # [] name of precipitation var
    "SpecHum_Field": "Q",                # [] name of specific humidity field
    "Evaporation_Field": "ET",           # [] name of evaporation field
    "PotentialEvaporation_Field": "E0",  # [] name of potential evap (for budyko evap calc)
    "BudykoOmega_Field": "omega_budyko", # [] name of budyko free parameter field

    # --- TOPOGRAPHY --------
    "orog_partition__TF": True,            # [True, False] use orog vs background precip components
    "topo_varname": "topo_m",              # [] name of topography variable
    "topo_mask_varname": "orog_mask",      # [] name of orography mask
    "slope_threshold_m_cell-1": 400,       # [m deg-1] slope threshold for orog mask
    "elev_threshold_m": 300,               # [m] elevation threshold for orog mask
    "UseOcean_in_OrogZonalMean__TF": True, # [True, False] do we use ocean precip in the zonal mean baseline



    # FORCING -------------------------------------------- # 
    # ---------------------------------------------------- #  

    # --- INSOLATION --------    
    "forcename_insol": "insolation",  # [] name of insolation forcing var
    "force_multiplier_insol": 1.,     # [] multiply insolation forcing by value
    "T_inertia_land": 0,              # [0-1] fraction of insol anomaly absorbed
    "T_inertia_ocean": 1,             # [0-1] fraction of insol anomaly absorbed
    # ... insol bounding box
    "force_bbox_insol": None,         # dict with "lat_bnds", "lon_bnds", and "cross_0_lon" bool

    # --- ALBEDO ------------
    "forcename_albedo": "albedo",           # [] name of albedo forcing var 
    "force_multiplier_albedo": 1.,          # [] multiply insolation forcing by value
    "TOAalbedo_sensitivity_factor": 0.7,    # [] fraction of surface albedo anom reaching TOA (Boos & Korty, 2016)

    # --- ARBITRARY ---------
    "forcename_arbitrary": "forcing",    # [] name of arbitrary forcing var
    "force_multiplier_arbitrary": 1.,    # [] multiply insolation forcing by value
    "arbitrary_forcing_Wm2": 0.,         # [W m-2] arbitrary forcing
    "arbitrary_forcing_landOnly": True,  # [True, False] apply forcing only over land?
    "force_bbox_arbitrary": None,        # dict with "lat_bnds", "lon_bnds", and "cross_0_lon" bool



    # ENERGY BALANCE MODEL ------------------------------- #
    # ---------------------------------------------------- #

    # --- EFE ---------------
    "EFE_absmax": 35,                    # [deg] max abslat for EFE
    "heaviside_deglim_V": 20,            # [deg] p/m latitude degrees around EFE where shifts are allowed
    "EFE_threshold": 0.1e8,              # [W m-1] V div max allowable for EFE
    "EFE_threshold_step": 0.02e8,        # [W m-1] increase EFE threshold by this much at a time if no EFE found
    "EFE_threshold_max": 0.5e8,          # [W m-1] max V div allowable for EFE *after* steps
    "DualEFE_ChiRangeForMean": 0.2e8,    # [W m-1] if two EFEs per some lon, take the mean if the Chis at EFE are within this range of each other, otherwise take the EFE @ lowest Chi

    # --- EFPM --------------
    "EFPM_nLinesMax_CTRL": 2,              # [int] maximum number of EFPM lines to define for control case
    "EFPM_nLinesMax_MATCH": 3,             # [int] max number of EFPM lines allowed for searching for a new match
    "EFPM_lat_absmax": 45,                 # [deg] max abslat for defining EFPM
    "EFPM_stitch_search_radius_deg": 17,   # [deg] abslon range to search for EFPM to stitch broken contours
    "EFPMs_Match_ProximityLimit": 50,      # [deg] lon range where EFPM in new clim state can be matched w/ EFPM from old clim state
    "EFPM_removeHorizontalLines_TF": True, # [True, False] whether to remove EFPM lines that surpass a horizontal angle threshold (likely EFE)
    "EFPM_HorizontalLine_angle": 30,       # [deg] remove lines shallower than this angle
    "EFPM_ChiThresholdFactor": 0.2,        # [0-1] fraction of minimum chi value for max allowable chi near EFPM
    "EFPM_MinLatSpan": 20,                 # [deg] EFPM must span this many degrees, latitude
    "heaviside_deglim_U": 70,              # [deg] p/m longitude degrees around EFPM where shifts are allowed
    "Lat_limit_of_zonal_P-shift": 45,      # [deg] max abslat for allowing precipitation shift to modify P



    # ISOTOPE ATTENUATION MODEL -------------------------- #
    # ---------------------------------------------------- #
    "SolveIsotopes": True,                  # [True, False] when false, only the EBM is run and all options in this section below have no effect

    # --- TAU-BAR MAP -------
    "Format_for_where_to_compute_tau-bar": "bbox",   # [] compute tau-bar following "bbox" or "coord_list"
    "CoordList_filename": "scratch_streamlines2save.csv",      # [] name of file w/ coords for computing iso. must be csv w/ lat, lon columns and in 'inputs' dir
    "Tau-bar_ComputeLandOnly": False,                       # [True, False] do not compute tau-bar over ocean?
    "Tau-bar_bbox_LatRange": [-36, 36],                     # [minlat, maxlat] limits for computing tau
    "Tau-bar_bbox_LonRange": [0, 360],                     # [minlon, maxlon] limits for computing tau
    "Tau-bar_bbox_Resolution": [2, 2],                 # [ndeg lat, ndeg lon] number of latitude and lnogitude degrees for sampling bbox (will get interpolated to cesm grid, but coarser grid saves compute time here)

    # --- ATTENUATION PARS --
    "Tau-bar_streamline-dx_km": 14,          # [km] interpolation step for transport vectors
    "Tau-bar_streamline-MaxTau": 8,          # [tau] maximum upwind tau before we deem attenuation negligible
    "Tau-bar_streamline-MaxDist_km": 25000,  # [km] maximum upwind distance for tau integration 
    "n_samples_TauBar_Land": 10,             # [int] number of steps to sample for getting E-source weighted tau-bar (for land evap correction)

    # --- TAU STREAMLINES ---
    "Collect_StreamlineDat_For-": "some",          # ["all", "some", "none"] how many points to collect streamline data for ("some" requires .csv with lat, lon, collect_streamline columns)
    "StreamlineSave_Coarsener": 10,                # [int] factor decrease in streamline resolution (*dx)

    # --- TAU DECOMPOSITION -
    # by space (!runs slower if we don't save streamlines where tau is computed... careful not to overload memory by saving too many though!)
    "Tau-bar_decomp_local_v_regional_LocalEvapEffect": False,  # [True, False] whether to do local v regional output with local P effect on local moisture change
    "Tau-bar_decomp_local_v_regional_UpwindEffect": True,      # [True, False] whether to do local v regional output with local P impacting remote moisture
    "Tau-bar_localThreshold": 1000,                  # [km] distance from sink point where we consider climatology to be local
    "Tau-bar_solve_halfEffect_dist_localEvap": False,    # [True, False] estimate the distance from sink within which half the Dtaubar change occurs when only local evap experiences local P
    "Tau-bar_solve_fracEffect_dist_streamline": True,    # [True, False] estimate the distance from sink within which half the Dtaubar change occurs when remote evap is allowed to be attenuated by local P change
    # by climatology (!will run slower!)
    "Tau-bar_decomp_E_L_s": False,                   # [True, False] whether to decompose tau-bar into changes due to Evap, Length scale, and flow field (requires 3x more streamline computations per timestep, so this slows things down a lot!)
    # interannual (True) or seasonal (False)
    "ELs_initstate_SameYrSlice": True,              # [True, False] when defining the baseline (init) state, do hold yr_slice constant and compare across time (i.e. time=1 mjjas compared to all other mjjas), or do we compare yr_slice within each tstep (i.e. at some tstep, mjjas compared to ndjfm)
    # fraction of Dtaubar to explain by "local" change
    "Tau-bar_DtauFraction": 0.75,          # [value 0-1] fraction of total Dtaubar that we want to attribute to "local" change, computing the distance before this fraction is achieved
    # Save composite streamlines too?
    "Tau-bar_save_local_v_regional_streamlines": True,  # [True, False] whether we store the local-only and regional-only streamline data 

    # CONSTANTS ------------------------------------------ # 
    # ---------------------------------------------------- #
    "Lv": 2.54e6,                # [J kg-1] latent heat of water vaporization (ignoring slight T dependency)
    "g": 9.81,                   # [m s-2] grav. acceleration



    # TIMING --------------------------------------------- # 
    # ---------------------------------------------------- #
    # ... these are defined once, then updated throughout 
    #     the run for print-out purposes
    "ts_year": 0,                  # [yr bp] timestep year
    "slice_of_year": '',           # [] ** will be updated from 'slices_to_solve' above

}   
# %%
