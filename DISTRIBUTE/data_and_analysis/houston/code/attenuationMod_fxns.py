# %% 
# ---------------------------------------------- # 
# iEBM FUNCTIONS                                 # 
# ---------------------------------------------- # 

# as v1, but with functionality for seasons 

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits
import glob
import pandas as pd
import xarray as xr
import cartopy
from cartopy.util import add_cyclic_point
import cartopy.crs as ccrs
import scipy.stats
import scipy.integrate
from scipy.interpolate import griddata  # for re-gridding np arrays in tau-bar fxn
import os
# import windspharm
import math
# import spharm
# import _spherepack
import dask
import warnings
# from windspharm.standard import VectorWind
from scipy.ndimage import distance_transform_edt # for getting distance to efe / efpm
from skimage import measure  # for contours
from scipy.spatial import distance # for efpm line stitching
from scipy.interpolate import RegularGridInterpolator  # for tau-bar computation
import warnings
import random

import accessory_fxns as fxn

# %% 


# ----------------------------------------------------------------
class Hydroclim:

    # [1] Partition orographic and non-orographic components
    def orog_partition(dataset, apply_to_var="PRECT",  **kwargs):
        # ... check if dataset already has topography ------------
        topo_var = kwargs.get("topo_varname")
        # provide name if none 
        if topo_var is None:
            topo_var = 'topo_m'
        # search for topo_var in dataset
        topo_var_check = topo_var in dataset.variables
        # check if it has a topo mask
        topo_mask = kwargs.get("topo_mask_varname")
        # provide name if none
        if topo_mask is None:
            topo_mask = 'orog_mask'
        # search for topo_var in dataset
        topo_mask_check = topo_mask in dataset.variables
        # --------------------------------------------------------

        # --- DEFINE SLOPE AND ELEV THRESHOLDS
        slope_threshold = kwargs.get("slope_threshold_m_cell-1")
        elev_threshold = kwargs.get("elev_threshold_m")
        if slope_threshold is None:
            slope_threshold = 400
        if elev_threshold is None:
            elev_threshold = 300

        # read in other topo file if needed
        if not topo_var_check:
            topo_file = kwargs.get("topo_file")
            if not topo_file:
                raise ValueError("I couldn't find topography in the input dataset, nor a specified path/file!")
            dat = xr.open_dataset(topo_file)
            # remove time dimension if exists
            if "time" in dat.coords:
                dat = dat.mean('time')

            # --- Extract terrestrial elev
            land_only = np.where(dat['data'] >= 0, dat['data'], 0)
            dat[topo_var] = xr.DataArray(land_only, dims=('lat', 'lon'),
                                            coords={'lat': dat.lat.values, 'lon': dat.lon.values},
                                            name=topo_var)
            
            # Interpolate and fill NAs
            dataset[topo_var] = dat[topo_var].interp_like(dataset, method = 'linear')
            # Fill extrapolated points with nearest value
            dataset[topo_var] = dataset[topo_var].ffill(dim='lon').bfill(dim='lon')
            
        # compute the orography mask if needed
        if not topo_mask_check:
            
            # ! Note ! Could add more methods here in future...

            # --- GET GRADIENT
            # Compute the slope using numpy.gradient
            elevation = dataset[topo_var].values
            lat = dataset['lat'].values
            lon = dataset['lon'].values
            slope_y, slope_x = np.gradient(elevation, lat, lon)

            # Combine the slopes into a single array
            slope = np.sqrt(slope_y ** 2 + slope_x ** 2)

            # Add the slope data to the xarray dataset
            dataset['slope'] = xr.DataArray(slope, dims=('lat', 'lon'), name='slope')

            
            mask = np.where(np.logical_and(dataset['slope'] >= slope_threshold, dataset[topo_var] > elev_threshold), 
                                  1, 0)
            dataset[topo_mask] = xr.DataArray(mask, dims=('lat', 'lon'),
                                                coords={'lat': dataset.lat.values, 'lon': dataset.lon.values},
                                                name=topo_mask)

        # --- DEFINE BACKGROUND VS MOVABLE NAMES
        if apply_to_var == "PRECT":
            var_name = kwargs.get("Inputfile_name_Precip")
            if var_name is None:
                var_name = apply_to_var
        
        # --- DEFINE BACKGROUND VS MOVABLE NAMES
        elif apply_to_var == "Q":
            var_name = kwargs.get("Inputfile_name_SpecHum")
            if var_name is None:
                var_name = apply_to_var

        # --- DEFINE BACKGROUND VS MOVABLE NAMES
        elif apply_to_var == "UQ":
            var_name = kwargs.get("Inputfile_name_ZonalQFlux")
            if var_name is None:
                var_name = apply_to_var
        
        # --- DEFINE BACKGROUND VS MOVABLE NAMES
        elif apply_to_var == "VQ":
            var_name = kwargs.get("Inputfile_name_MeridQFlux")
            if var_name is None:
                var_name = apply_to_var

        else:
            var_name = apply_to_var
        
        var_name_fixed = var_name + "_fixed"
        var_name_movable = var_name + "_movable"
        
        
        
        # do we include ocean in zonal mean?
        ocean_plus_land = kwargs.get("UseOcean_in_OrogZonalMean__TF")
        if ocean_plus_land:
            zonal_mean_var = dataset[var_name].mean(dim='lon')  
        else:
            zonal_mean_var = dataset[var_name].where(dataset['LANDFRAC']).mean(dim='lon')  
        
        # compute anomaly to zonal mean
        anomaly_var = dataset[var_name] - zonal_mean_var
        
        # store in dataset and return
        v_fixed = np.where(dataset[topo_mask], anomaly_var, 0)
        dataset[var_name_fixed] = xr.DataArray(v_fixed, dims=('lat', 'lon'),
                                                coords={'lat': dataset.lat.values, 'lon': dataset.lon.values},
                                                name=var_name_fixed)
        dataset[var_name_movable] = dataset[var_name] - dataset[var_name_fixed]

        return dataset
    

    

    
    






# ----------------------------------------------------------------
class Isotopes:

    # [1] COMPUTE TAU BAR -- no decomposing the components
    @staticmethod
    def taubar_noDecomposition(ds, iso_coords=None, **kwargs):
        # skip if we're in no EBM mode
        if (kwargs.get('runtype') is not None) & ("noIsotopes" in kwargs.get('runtype')):
            return ds, None
        # -----------------------------------------------
        # FIGURE OUT WHICH CELLS TO COMPUTE -----------------------
        # ... should be 'bbox' or 'coord_list'
        compute_method = kwargs.get('Format_for_where_to_compute_tau-bar')
        streamline_collect = kwargs.get("Collect_StreamlineDat_For-")
        land_only = kwargs.get('Tau-bar_ComputeLandOnly')
        n_samples = kwargs.get('n_samples_TauBar_Land')
        
        if land_only is None:
            land_only = True
        
        if n_samples is None:   # number of land steps to sample for evap correction of taubar
            n_samples = 5

        if streamline_collect is None:
            streamline_collect = "none"  # don't output streamlines
        if compute_method is None:
            compute_method = 'bbox'
        
        if compute_method == 'bbox':
            lat_range = kwargs.get('Tau-bar_bbox_LatRange')
            lon_range = kwargs.get('Tau-bar_bbox_LonRange')
            res = kwargs.get('Tau-bar_bbox_Resolution')
            if lat_range is None:
                lat_range = [-35, 35]
            if lon_range is None:
                lon_range = [0, 360]
            if res is None:
                res_adj = False
            else:
                res_adj = True
        
        n_samples_init = n_samples # we may have to update later
        # CHECK FOR REQUIRED VARS ---------------------------------
        vars_required = ['UQ', 'VQ', kwargs.get("Evaporation_Field"), 
                         kwargs.get("Precipitation_Field"), "Fmag"]
        fxn.var_check(vars_required, ds)

        # SET TAU GRID AND TRACKER --------------------------------
        year_printout = kwargs.get('ts_year')
        slice_printout = kwargs.get('slice_of_year')
        yrslc_remain_printout = kwargs.get('yr_slices_remaining')
        if year_printout is None:
            year_printout = "unknown..."
        if slice_printout is None:
            slice_printout = "unknown..."
        if yrslc_remain_printout is None:
            yrslc_remain_printout = "unknown..."

        dx = kwargs.get("Tau-bar_streamline-dx_km")         # represents the distance increment along the streamline for calculating tau_bar (km)
        if dx is None:
            dx = 8
        Dx = dx/111 # the increment in degrees latitude (and, at the equator, longitude).
        taumax = kwargs.get("Tau-bar_streamline-MaxTau") # maximum tau to integrate along upstream path
        if taumax is None:
            taumax = 8
        dmax = kwargs.get("Tau-bar_streamline-MaxDist_km")
        if dmax is None:
            dmax = 25000         # maximum distance to integrate along the upstream path (km)
        Nmax = math.ceil(dmax/dx) # maximum number of calculations along the streamline
        # ---------------------------------------------------------

        # --- SET UP LOOP
        # indices
        lat = ds.lat
        lon = ds.lon

        # initialize arrays
        P, E, UQ, VQ, Fmag, lfrac = fxn.tau_arrayInitialize(ds, lat, lon, **kwargs)
        # get interpolation functions
        Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = fxn.tau_interpolatorInitialize(lon, lat, E, P, UQ, VQ, lfrac)
        
        # set bbox grid
        if compute_method == 'bbox':
            if res_adj:
                lat_temp = np.arange(-90, 90+1, 1, dtype=int)
                lon_temp = np.arange(0, 360+1, 1, dtype=int)
                lat = lat_temp[::res[0]]
                lon = lon_temp[::res[1]]
                # Define the latitude/longitude range you want to select
                lat_bool = np.logical_and(lat >= lat_range[0], lat <= lat_range[1])
                lon_bool = np.logical_and(lon >= lon_range[0], lon <= lon_range[1])
            else:
                # Define the latitude/longitude range you want to select
                lat_bool = np.logical_and(lat.lat.values >= lat_range[0], lat.lat.values <= lat_range[1])
                lon_bool = np.logical_and(lon.lon.values >= lon_range[0], lon.lon.values <= lon_range[1])


        # empty tau-bar array
        tau_bar = np.full((lat.shape[0],lon.shape[0]), np.nan) # create empty output array
        tau_bar_Ewtd_land = np.full_like(tau_bar, np.nan) 
        moisture_dist_inland = np.full_like(tau_bar, np.nan) 
        land_frac_of_streamline = np.full_like(tau_bar, np.nan) 
        land_frac_of_evapSource = np.full_like(tau_bar, np.nan) 

        # -- LOOP THROUGH LATITUDES -- 
        # BBOX LOOP --------------------------------------------------------------------
        if compute_method == 'bbox': 
            # find where we save streamline (if at all)
            if streamline_collect == "some":
                iso_coords_collect = iso_coords.loc[iso_coords['collect_streamline']=="Y",]
                # get in 0-360
                iso_coords_collect['lon2'] = np.where(iso_coords_collect['lon'] < 0, iso_coords_collect['lon'] % 360, iso_coords_collect['lon'])
                # get x and y indices of nearest lat/lon in grid (for storing result)
                if res_adj:
                    lat_dx_collect = np.argmin(np.abs(lat[:, np.newaxis] - np.array(iso_coords_collect['lat'])), axis=0)
                    lon_dx_collect = np.argmin(np.abs(lon[:, np.newaxis] - np.array(iso_coords_collect['lon2'])), axis=0)
                    latval = lat[lat_dx_collect]
                    lonval = lon[lon_dx_collect]
                    xsave = np.zeros(len(lonval), dtype=int); ysave = np.zeros_like(xsave)
                    for iso_row in range(len(iso_coords_collect)):
                        xsave[iso_row] = np.where(lon == lonval[iso_row])[0][0]
                        ysave[iso_row] = np.where(lat == latval[iso_row])[0][0]
                else:
                    latval = ds.sel(lat=iso_coords_collect['lat'], method='nearest').lat.values
                    lonval = ds.sel(lon=iso_coords_collect['lon2'], method='nearest').lon.values
                    xsave = np.zeros(len(lonval), dtype=int); ysave = np.zeros_like(xsave)
                    for iso_row in range(len(iso_coords_collect)):
                        xsave[iso_row] = np.where(ds.lon.values == lonval[iso_row])[0][0]
                        ysave[iso_row] = np.where(ds.lat.values == latval[iso_row])[0][0]

            # Get the indices 
            A_slice = np.where(lat_bool)[0]
            B_slice = np.where(lon_bool)[0]
            # initialize tracker
            idx_counter = 1
            n_cells = len(A_slice) * len(B_slice)
            for y in A_slice:
                # -- LOOP THROUGH LONGITUDES -- 
                for x in B_slice:
                    # coordinates at point of precipitation
                    lat0 = lat[y]
                    lon0 = lon[x]
                    # save sink value (above values are overwritten)
                    lat0_init = lat0
                    lon0_init = lon0
                    
                    # check whether to save streamline
                    if streamline_collect == "some":
                        # check lon, then lat
                        x_check = np.where(xsave == x)[0]
                        if x_check.size == 0:
                            save_streamline = False
                        else: # if we match an x, check if we match a y
                            if y in ysave[x_check]:
                                save_streamline = True
                            else:
                                save_streamline = False
                    elif streamline_collect == "all":
                        save_streamline = True

                    # check if it's on land
                    land_check = ds['LANDFRAC'].sel(lat=lat0, lon=lon0, method='nearest').values
                    if land_only and land_check == 0.:
                        moisture_dist_inland[y,x] = 0
                        # track progress -------------------------------------------------------------
                        complete_percent = round((idx_counter / n_cells)*100,2)
                        l1 = "...skipping tau over the ocean " + " -- " + str(complete_percent) + "%" + " complete"
                        l2 = "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ "
                        print(l1, l2, sep=os.linesep)
                        # ----------------------------------------------------------------------------
                        idx_counter += 1 # move to next index
                        continue
                    elif land_check < 0.6:
                        moisture_dist_inland[y,x] = 0

                    # compute tau and climatology across streamline
                    coast_step, nanindex, dist, dist_inland, tau, mu, wp, E0, P0, vq0, uq0, Fmag0, lfrac0, PminE, lat_save, lon_save = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)
                    # compute E source fraction terrestrial
                    tau_wtd_mean, land_Esource_frac, E_terrFrac = fxn.terrestrial_E_frac(land_check, wp, coast_step, tau, E0, nanindex, n_samples)
                    
                    
                    # reset n_samples
                    n_samples = n_samples_init  # probably unnecessary now that this is changed in a separate fxn, but can't hurt
                    # store wtd result
                    tau_bar_Ewtd_land[y,x] = E_terrFrac * tau_wtd_mean
                    # moisture dist inland
                    moisture_dist_inland[y,x] = dist_inland
                    # fraction of streamline over land 
                    land_frac_of_streamline[y,x] = dist_inland / np.nanmax(dist)
                    # fraction of evap source coming from land ET 
                    land_frac_of_evapSource[y,x] = land_Esource_frac

                    # construct the output dataframe ------------------
                    no_streamlines_yet = 'df_streamlines' not in locals()
                    if no_streamlines_yet:
                        df_streamlines, tau_x, E0_x, dist_x, mu_x = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                            save_streamline, dist, lat_save, lon_save, E0, Fmag0, P0, PminE, 
                                                                                            tau, mu, wp, lfrac0, nanindex, streamline_type='full_model', **kwargs)
                    else:
                        df_streamlines, tau_x, E0_x, dist_x, mu_x = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                        save_streamline, dist, lat_save, lon_save, E0, Fmag0, P0, PminE, 
                                                                                        tau, mu, wp, lfrac0, nanindex, df_streamlines, streamline_type='full_model', **kwargs)
                    # --------------------------------------------------

                    # DIAGNOSTIC - Turned off even if these lines are uncommmented (would need to output stuff in tau_streamline_point fxn)
                    # if get_taubar_streamline:
                    #     idx = np.where(dist[nanindex] <= check_last_km)[0]
                    #     for step in idx:
                    #         tau_temp = tau[step:] - tau[step] # start at zero (subtract all upwind)
                    #         tau_bar_streamline[step] = np.trapz(tau_temp*E0[step:]*math.e**(-tau_temp))/np.trapz(E0[step:]*math.e**(-tau_temp))

                    # resample to tau grid (instead of space) ----------
                    # (this is necessary to decompose E and L)
                    mu_taugrid, wp_taugrid, E0_taugrid, tausteps = fxn.dist_to_tau_coords(tau, nanindex, mu, wp, E0)

                    # TAU_BAR COMPUTATION
                    tau_bar[y,x] = np.trapz(E0_taugrid * (1/mu_taugrid) * tausteps * math.e**(-tausteps)) / np.trapz(E0_taugrid * (1/mu_taugrid) * math.e**(-tausteps))
                
                    
                    # track progress -------------------------------------------------------------
                    complete_percent = round((idx_counter / n_cells)*100,2)
                    l1 = "Year: " + year_printout + ", Slice: " + slice_printout + " -- " + str(complete_percent) + "%" + " complete"
                    l2 = "(" + yrslc_remain_printout + " Year-Slices remaining)"
                    l3 = "--------------------------------------------"
                    print(l1, l2, l3, sep=os.linesep)
                    # ----------------------------------------------------------------------------
                    
                    # adjust counter
                    idx_counter += 1
                    # print("complete")
            

            # restore to original resolution
            if res_adj:
                old_Y, old_X = np.meshgrid(lon, lat)
                new_y = ds.lat.values
                new_x = ds.lon.values
                new_Y, new_X = np.meshgrid(new_x, new_y)
                # update arrays
                tau_bar = griddata((old_Y.flatten(), old_X.flatten()), tau_bar.flatten(), (new_Y, new_X), method='linear')
                tau_bar_Ewtd_land = griddata((old_Y.flatten(), old_X.flatten()), tau_bar_Ewtd_land.flatten(), (new_Y, new_X), method='linear')
                moisture_dist_inland = griddata((old_Y.flatten(), old_X.flatten()), moisture_dist_inland.flatten(), (new_Y, new_X), method='linear')
                land_frac_of_streamline = griddata((old_Y.flatten(), old_X.flatten()), land_frac_of_streamline.flatten(), (new_Y, new_X), method='linear')
                land_frac_of_evapSource = griddata((old_Y.flatten(), old_X.flatten()), land_frac_of_evapSource.flatten(), (new_Y, new_X), method='linear')
                

        # COORD-LIST LOOP -----------------------------------------------------------------
        if compute_method == 'coord_list':     
            idx_counter = 1
            n_cells = len(iso_coords)
            # -- LOOP THROUGH LONGITUDES -- 
            for coord in range(n_cells):
                # coordinates at point of precipitation
                lat0 = iso_coords['lat'].iloc[coord]
                lon0 = iso_coords['lon'].iloc[coord]

                # check whether to save streamline
                if streamline_collect == "some":
                    if iso_coords['collect_streamline'].iloc[coord] == "Y":
                        save_streamline = True
                elif streamline_collect == "all":
                    save_streamline = True

                # convert to 0-360 if necessary
                if lon0 < 0:
                    lon0 = lon0 % 360

                # get x and y indices of nearest lat/lon in grid (for storing result)
                latval = ds.sel(lat=lat0, method='nearest').lat.values
                lonval = ds.sel(lon=lon0, method='nearest').lon.values
                x = np.where(ds.lon.values == lonval)[0][0]
                y = np.where(ds.lat.values == latval)[0][0]

                # check if it's on land
                land_check = ds['LANDFRAC'].sel(lat=lat0, lon=lon0, method='nearest').values
                if land_only and land_check == 0.:
                    moisture_dist_inland[y,x] = 0
                    # track progress -------------------------------------------------------------
                    complete_percent = round((idx_counter / n_cells)*100,2)
                    l1 = "...skipping tau over the ocean " + " -- " + str(complete_percent) + "%" + " complete"
                    l2 = "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ "
                    print(l1, l2, sep=os.linesep)
                    # ----------------------------------------------------------------------------

                    idx_counter += 1 # move to next index
                    continue
                elif land_check < 0.6:
                    moisture_dist_inland[y,x] = 0

                # compute tau and climatology across streamline
                coast_step, nanindex, dist, dist_inland, tau, mu, wp, E0, P0, vq0, uq0, Fmag0, lfrac0, PminE, lat_save, lon_save = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)
                # compute E source fraction terrestrial
                tau_wtd_mean, land_Esource_frac, E_terrFrac = fxn.terrestrial_E_frac(land_check, wp, coast_step, tau, E0, nanindex, n_samples)
                
                # reset n_samples
                n_samples = n_samples_init
                # store wtd result
                tau_bar_Ewtd_land[y,x] = E_terrFrac * tau_wtd_mean
                # moisture dist inland
                moisture_dist_inland[y,x] = dist_inland
                # fraction of land the streamline experiences
                land_frac_of_streamline[y,x] = dist_inland / np.nanmax(dist)
                # fraction of evap source coming from land ET 
                land_frac_of_evapSource[y,x] = land_Esource_frac


                # construct the output dataframe ------------------
                no_streamlines_yet = 'df_streamlines' not in locals()
                if no_streamlines_yet:
                    df_streamlines, tau_x, E0_x, dist_x, mu_x = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                        save_streamline, dist, lat_save, lon_save, E0, Fmag0, P0, PminE, 
                                                                                        tau, mu, wp, lfrac0, nanindex, **kwargs)
                else:
                    df_streamlines, tau_x, E0_x, dist_x, mu_x = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                    save_streamline, dist, lat_save, lon_save, E0, Fmag0, P0, PminE, 
                                                                                    tau, mu, wp, lfrac0, nanindex, df_streamlines, **kwargs)
                # --------------------------------------------------

                # if get_taubar_streamline:
                #     idx = np.where(dist[nanindex] <= check_last_km)[0]
                #     for step in idx:
                #         tau_temp = tau[step:] - tau[step] # start at zero (subtract all upwind)
                #         tau_bar_streamline[step] = np.trapz(tau_temp*E0[step:]*math.e**(-tau_temp))/np.trapz(E0[step:]*math.e**(-tau_temp))
                
                # resample to tau grid (instead of space) ----------
                # (this is necessary to decompose E and L)
                mu_taugrid, wp_taugrid, E0_taugrid, tausteps = fxn.dist_to_tau_coords(tau, nanindex, mu, wp, E0)

                # TAU_BAR COMPUTATION
                tau_bar[y,x] = np.trapz(E0_taugrid * (1/mu_taugrid) * tausteps * math.e**(-tausteps)) / np.trapz(E0_taugrid * (1/mu_taugrid) * math.e**(-tausteps))



                # track progress -------------------------------------------------------------
                complete_percent = round((idx_counter / n_cells)*100,2)
                l1 = "Year: " + year_printout + ", Slice: " + slice_printout + " -- " + str(complete_percent) + "%" + " complete"
                l2 = "(" + yrslc_remain_printout + " Year-Slices remaining)"
                l3 = "--------------------------------------------"
                print(l1, l2, l3, sep=os.linesep)
                # ----------------------------------------------------------------------------

                
                # adjust counter
                idx_counter += 1
                
        

        # ADD RESULT TO DS
        tau_name = "tau_bar"
        moisture_dist_name = "moisture_dist_inland"
        tau_bar_EWL_name = "tau_bar_wtdLandEvap"
        landfrac_streamline_name = "streamline_frac_land"
        land_frac_of_evapSource_name = "Esource_frac_land"
        # get to data array
        ds[tau_name] = xr.DataArray(tau_bar[:,:], dims=('lat', 'lon'),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values},
                                    name=tau_name)
        ds[moisture_dist_name] = xr.DataArray(moisture_dist_inland[:,:], dims=('lat', 'lon'),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values},
                                    name=moisture_dist_name)
        ds[tau_bar_EWL_name] = xr.DataArray(tau_bar_Ewtd_land[:,:], dims=('lat', 'lon'),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values},
                                    name=tau_bar_EWL_name)
        ds[landfrac_streamline_name] = xr.DataArray(land_frac_of_streamline[:,:], dims=('lat', 'lon'),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values},
                                    name=landfrac_streamline_name)
        ds[land_frac_of_evapSource_name] = xr.DataArray(land_frac_of_evapSource[:,:], dims=('lat', 'lon'),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values},
                                    name=land_frac_of_evapSource_name)

        # return result
        return ds, df_streamlines



    # [2] COMPUTE TAU BAR -- with decomposition
    @staticmethod
    def taubar_withDecomposition(ds, iso_coords=None, **kwargs):
        # skip if we're in no EBM mode
        if (kwargs.get('runtype') is not None) & ("noIsotopes" in kwargs.get('runtype')):
            return ds, None
        
        # -----------------------------------------------
        # SET UP THE DECOMPOSITION
        baseline = kwargs.get('Tau-bar_decomposition_baseline')
        if kwargs.get('ELs_initstate_SameYrSlice'):
            decomp_dim = "time"
            # if there's a yr_slice dim, pull it out
            if 'yr_slice' in ds.dims:
                if len(ds.yr_slice) > 1:
                    raise ValueError("Too many yr_slices in dataset to decompose taubar over time. Must have 1 or no yr_slices. Maybe 'ELs_initstate_SameYrSlice' in initialize.py should be False?")
                else:
                    ds = ds.isel(yr_slice = 0)   # just take the first (only) one
            # collect or assign baseline if needed
            if (baseline is None) or (baseline not in ds[decomp_dim]): 
                baseline = np.min(ds.time.values)

        elif not kwargs.get('ELs_initstate_SameYrSlice'):
            decomp_dim = "yr_slice"
            # if there's a time dim, pull it out
            if 'time' in ds.dims:
                if len(ds.time) > 1:
                    raise ValueError("Too many time slices in dataset to decompose taubar over yr_slice. Must have 1 or no time slices. Maybe 'ELs_initstate_SameYrSlice' in initialize.py should be True?")
                else:
                    ds = ds.isel(time = 0)   # just take the first (only) one
            # collect or assign baseline if needed
            if (baseline is None) or (baseline not in ds[decomp_dim]): 
                baseline = ds.yr_slice.values[0]
        # index for baseline case
        baseline_idx = np.where(ds[decomp_dim].values == baseline)[0][0]


        # FIGURE OUT WHICH CELLS TO COMPUTE -----------------------
        # ... should be 'bbox' or 'coord_list'
        compute_method = kwargs.get('Format_for_where_to_compute_tau-bar')
        streamline_collect = kwargs.get("Collect_StreamlineDat_For-")
        land_only = kwargs.get('Tau-bar_ComputeLandOnly')
        n_samples = kwargs.get('n_samples_TauBar_Land')
        localEvapEff = kwargs.get('Tau-bar_decomp_local_v_regional_LocalEvapEffect')
        UpwindEff = kwargs.get('Tau-bar_decomp_local_v_regional_UpwindEffect')
        save_upwindEff_streamline = kwargs.get('Tau-bar_save_local_v_regional_streamlines')
        local_threshold = kwargs.get('Tau-bar_localThreshold')
        tau_ELs_decompose = kwargs.get('Tau-bar_decomp_E_L_s')
        Dtau_50_solve_localEvap = kwargs.get('Tau-bar_solve_fracEffect_dist_localEvap')
        Dtau_50_solve_strmline = kwargs.get('Tau-bar_solve_fracEffect_dist_streamline')
        Dtau_frac_upwindSolve = kwargs.get('Tau-bar_DtauFraction')

        if land_only is None:
            land_only = True
        
        if n_samples is None:   # number of land steps to sample for evap correction of taubar
            n_samples = 5

        if streamline_collect is None:
            streamline_collect = "none"  # don't output streamlines
        elif local_threshold is None:
            local_threshold = 750   # doesn't do anything if tau_local_v_regional_* is False
        if tau_ELs_decompose is None:
            tau_ELs_decompose = False
        if localEvapEff is None:
            localEvapEff = False
        if UpwindEff is None:
            UpwindEff = False
        if Dtau_50_solve_localEvap is None:
            Dtau_50_solve_localEvap = False
        if Dtau_50_solve_strmline is None:
            Dtau_50_solve_strmline = False
        if Dtau_frac_upwindSolve is None:
            Dtau_frac_upwindSolve = 1-(1/math.e)
        if compute_method is None:
            compute_method = 'bbox'
        
        if compute_method == 'bbox':
            lat_range = kwargs.get('Tau-bar_bbox_LatRange')
            lon_range = kwargs.get('Tau-bar_bbox_LonRange')
            res = kwargs.get('Tau-bar_bbox_Resolution')
            if lat_range is None:
                lat_range = [-35, 35]
            if lon_range is None:
                lon_range = [0, 360]
            if res is None:
                res_adj = False
            else:
                res_adj = True
        
        n_samples_init = n_samples # we may have to update later
        # CHECK FOR REQUIRED VARS ---------------------------------
        vars_required = ['UQ', 'VQ', kwargs.get("Evaporation_Field"), 
                            kwargs.get("Precipitation_Field"), "Fmag"]
        fxn.var_check(vars_required, ds)

        # SET TAU GRID AND TRACKER --------------------------------
        this_runName = kwargs.get('run_name')
        year_printout = kwargs.get('ts_year')
        slice_printout = kwargs.get('slice_of_year')
        yrslc_remain_printout = kwargs.get('yr_slices_remaining')
        total_steps = kwargs.get('total_steps')
        this_step = kwargs.get('this_step')
        if year_printout is None:
            year_printout = "unknown..."
        if slice_printout is None:
            slice_printout = "unknown..."
        if yrslc_remain_printout is None:
            yrslc_remain_printout = "unknown..."

        dx = kwargs.get("Tau-bar_streamline-dx_km")         # represents the distance increment along the streamline for calculating tau_bar (km)
        if dx is None:
            dx = 8
        Dx = dx/111 # the increment in degrees latitude (and, at the equator, longitude).
        taumax = kwargs.get("Tau-bar_streamline-MaxTau") # maximum tau to integrate along upstream path
        if taumax is None:
            taumax = 8
        dmax = kwargs.get("Tau-bar_streamline-MaxDist_km")
        if dmax is None:
            dmax = 25000         # maximum distance to integrate along the upstream path (km)
        Nmax = math.ceil(dmax/dx) # maximum number of calculations along the streamline
        # ---------------------------------------------------------

        # --- SET UP LOOP
        # indices
        lat = ds.lat
        lon = ds.lon
        kwargs['decomp'] = ds[decomp_dim]

        # initialize arrays
        P, E, UQ, VQ, Fmag, lfrac = fxn.tau_arrayInitialize(ds, lat, lon, **kwargs)
        
        # get interpolation functions for baseline case 
        Efiti, Pfiti, uqfiti, vqfiti, Fmagfiti, lfracfiti = fxn.tau_interpolatorInitialize(lon, lat, E[:,:,baseline_idx], P[:,:,baseline_idx], 
                                                                                    UQ[:,:,baseline_idx], VQ[:,:,baseline_idx], lfrac[:,:,baseline_idx])
        
        # handle resolution if needed
        if compute_method == 'bbox':
            if res_adj:
                lat_temp = np.arange(-90, 90+1, 1, dtype=int)
                lon_temp = np.arange(0, 360+1, 1, dtype=int)
                lat = lat_temp[::res[0]]
                lon = lon_temp[::res[1]]
                # Define the latitude/longitude range you want to select
                lat_bool = np.logical_and(lat >= lat_range[0], lat <= lat_range[1])
                lon_bool = np.logical_and(lon >= lon_range[0], lon <= lon_range[1])
            else:
                # Define the latitude/longitude range you want to select
                lat_bool = np.logical_and(lat.lat.values >= lat_range[0], lat.lat.values <= lat_range[1])
                lon_bool = np.logical_and(lon.lon.values >= lon_range[0], lon.lon.values <= lon_range[1])
        


        # empty tau-bar array
        tau_bar = np.full((lat.shape[0],lon.shape[0], len(kwargs.get('decomp'))), np.nan) # create empty output array
        tau_bar_Ewtd_land = np.full_like(tau_bar, np.nan) 
        moisture_dist_inland = np.full_like(tau_bar, np.nan) 
        land_frac_of_streamline = np.full_like(tau_bar, np.nan) 
        land_frac_of_evapSource = np.full_like(tau_bar, np.nan) 
        if localEvapEff:
            Dtau_bar_regionalEff_localEvap = np.full_like(tau_bar, np.nan)
            Dtau_bar_localEff_localEvap = np.full_like(tau_bar, np.nan)
        if UpwindEff:
            Dtau_bar_regionalEff_upwind = np.full_like(tau_bar, np.nan)
            Dtau_bar_localEff_upwind = np.full_like(tau_bar, np.nan)
        if tau_ELs_decompose:
            Dtau_bar_newE = np.full_like(tau_bar, np.nan)
            Dtau_bar_newL = np.full_like(tau_bar, np.nan)
            Dtau_bar_newPath = np.full_like(tau_bar, np.nan)
        if localEvapEff or UpwindEff or tau_ELs_decompose:  # with *any* decomposition, save the Dtau_bar total for convenience
            Dtau_bar = np.full_like(tau_bar, np.nan)
        if Dtau_50_solve_localEvap:
            Dtau_bar_50_dist_localEvap = np.full_like(tau_bar, np.nan)
        if Dtau_50_solve_strmline:
            Dtau_bar_50_dist_strmline = np.full_like(tau_bar, np.nan)


        # -- LOOP THROUGH LATITUDES -- 
        # BBOX LOOP --------------------------------------------------------------------
        if compute_method == 'bbox': 
            # find where we save streamline (if at all)
            if streamline_collect == "some":
                iso_coords_collect = iso_coords.loc[iso_coords['collect_streamline']=="Y",]
                # get in 0-360
                iso_coords_collect['lon2'] = np.where(iso_coords_collect['lon'] < 0, iso_coords_collect['lon'] % 360, iso_coords_collect['lon'])
                # get x and y indices of nearest lat/lon in grid (for storing result)
                if res_adj:
                    lat_dx_collect = np.argmin(np.abs(lat[:, np.newaxis] - np.array(iso_coords_collect['lat'])), axis=0)
                    lon_dx_collect = np.argmin(np.abs(lon[:, np.newaxis] - np.array(iso_coords_collect['lon2'])), axis=0)
                    latval = lat[lat_dx_collect]
                    lonval = lon[lon_dx_collect]
                    xsave = np.zeros(len(lonval), dtype=int); ysave = np.zeros_like(xsave)
                    for iso_row in range(len(iso_coords_collect)):
                        xsave[iso_row] = np.where(lon == lonval[iso_row])[0][0]
                        ysave[iso_row] = np.where(lat == latval[iso_row])[0][0]
                else:
                    latval = ds.sel(lat=iso_coords_collect['lat'], method='nearest').lat.values
                    lonval = ds.sel(lon=iso_coords_collect['lon2'], method='nearest').lon.values
                    xsave = np.zeros(len(lonval), dtype=int); ysave = np.zeros_like(xsave)
                    for iso_row in range(len(iso_coords_collect)):
                        xsave[iso_row] = np.where(ds.lon.values == lonval[iso_row])[0][0]
                        ysave[iso_row] = np.where(ds.lat.values == latval[iso_row])[0][0]

            # Get the indices 
            A_slice = np.where(lat_bool)[0]
            B_slice = np.where(lon_bool)[0]
            decomp_slice = range(len(ds[decomp_dim])-1)  # subtract 1 because we compute base case outside the decomposition dim loop
            case_idxs = np.where(ds[decomp_dim].values != baseline)[0]   # all indices except baseline
            # initialize tracker
            idx_counter = 1
            readout_counter = 1
            n_cells = len(A_slice) * len(B_slice) * len(decomp_slice)
            for y in A_slice:
                # -- LOOP THROUGH LONGITUDES -- 
                for x in B_slice:
                    # coordinates at point of precipitation
                    lat0 = lat[y]
                    lon0 = lon[x]
                    # save sink value (above values are overwritten)
                    lat0_init = lat0
                    lon0_init = lon0
                    
                    # check whether to save streamline
                    if streamline_collect == "some":
                        # check lon, then lat
                        x_check = np.where(xsave == x)[0]
                        if x_check.size == 0:
                            save_streamline = False
                        else: # if we match an x, check if we match a y
                            if y in ysave[x_check]:
                                save_streamline = True
                            else:
                                save_streamline = False
                    elif streamline_collect == "all":
                        save_streamline = True

                    # check if it's on land (use first value b/c should all be the same (over time or yr_slice dim))
                    land_check = ds['LANDFRAC'].sel(lat=lat0, lon=lon0, method='nearest').values[0]
                    if land_only and land_check == 0.:
                        moisture_dist_inland[y,x,:] = 0
                        # track progress -------------------------------------------------------------
                        complete_percent = round((readout_counter / n_cells)*100,2)
                        l1 = "...skipping tau over the ocean " + " -- " + str(complete_percent) + "%" + " complete"
                        l2 = "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ "
                        print(l1, l2, sep=os.linesep)
                        # ----------------------------------------------------------------------------
                        readout_counter += 1 # move to next index
                        continue
                    elif land_check < 0.6:
                        moisture_dist_inland[y,x,baseline_idx] = 0

                    # compute tau and climatology across streamline
                    coast_stepi, nanindexi, disti, dist_inlandi, taui, mui, wpi, E0i, P0i, vq0i, uq0i, Fmag0i, lfrac0i, PminEi, lat_savei, lon_savei = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                Efiti, Pfiti, vqfiti, uqfiti, Fmagfiti, lfracfiti)
                    # compute E source fraction terrestrial
                    tau_wtd_meani, land_Esource_fraci, E_terrFraci = fxn.terrestrial_E_frac(land_check, wpi, coast_stepi, taui, E0i, nanindexi, n_samples)
                    
                    # construct the output dataframe ------------------
                    # set "df_streamlines" to none each time so it just gives the current streamline and we add decomp_dim and concatenate after the fact
                    no_streamlines_yet = True   # hold always true for this version
                    df_streamlinesi, _, _, _, _ = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                        save_streamline, disti, lat_savei, lon_savei, E0i, Fmag0i, P0i, PminEi, 
                                                                                        taui, mui, wpi, lfrac0i, nanindexi, streamline_type='full_model', df_streamlines=None, **kwargs)
                    df_streamlinesi[decomp_dim] = baseline
                    if 'df_streamlines' not in locals():
                        df_streamlines = df_streamlinesi.copy()
                    else:
                        df_streamlines = pd.concat([df_streamlines, df_streamlinesi], ignore_index=True)
                    # --------------------------------------------------

                    # reset n_samples
                    n_samples = n_samples_init  # probably unnecessary now that this is changed in a separate fxn, but can't hurt
                    # store wtd result
                    tau_bar_Ewtd_land[y,x,baseline_idx] = E_terrFraci * tau_wtd_meani
                    # moisture dist inland
                    moisture_dist_inland[y,x,baseline_idx] = dist_inlandi
                    # fraction of streamline over land 
                    land_frac_of_streamline[y,x,baseline_idx] = dist_inlandi / np.nanmax(disti)
                    # fraction of evap source coming from land ET 
                    land_frac_of_evapSource[y,x,baseline_idx] = land_Esource_fraci

                    # resample to tau grid (instead of space) ----------
                    # (this is necessary to decompose E and L)
                    mu_taugridi, wp_taugridi, E0_taugridi, taustepsi = fxn.dist_to_tau_coords(taui, nanindexi, mui, wpi, E0i)

                    # SOLVE TAU-BAR
                    tau_bar[y,x,baseline_idx] = np.trapz(E0_taugridi * (1/mu_taugridi) * taustepsi * math.e**(-taustepsi)) / np.trapz(E0_taugridi * (1/mu_taugridi) * math.e**(-taustepsi))
                    taubar_init = tau_bar[y,x,baseline_idx] 


                    # -- LOOP THROUGH DECOMPOSITION DIM -- 
                    for ddim in decomp_slice:
                        thisidx = case_idxs[ddim]
                        thisval = ds[decomp_dim].values[thisidx]
                        # build the grid interpolator 
                        # get interpolation functions for baseline case 
                        Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = fxn.tau_interpolatorInitialize(ds.lon.values, ds.lat.values, E[:,:,thisidx], P[:,:,thisidx], 
                                                                                                    UQ[:,:,thisidx], VQ[:,:,thisidx], lfrac[:,:,thisidx])
                        # compute tau and climatology across streamline
                        coast_step, nanindex, dist, dist_inland, tau, mu, wp, E0, P0, vq0, uq0, Fmag0, lfrac0, PminE, lat_save, lon_save = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                    Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)
                        # compute E source fraction terrestrial
                        tau_wtd_mean, land_Esource_frac, E_terrFrac = fxn.terrestrial_E_frac(land_check, wp, coast_step, tau, E0, nanindex, n_samples)
                        
                        # update streamlines
                        # construct the output dataframe ------------------
                        # set "df_streamlines" to none each time so it just gives the current streamline and we add decomp_dim and concatenate after the fact
                        no_streamlines_yet = True   # hold always true for this version
                        tmp_streamlines, _, _, _, _ = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                            save_streamline, dist, lat_save, lon_save, E0, Fmag0, P0, PminE, 
                                                                                            tau, mu, wp, lfrac0, nanindex, df_streamlines=None, streamline_type='full_model', **kwargs)
                        tmp_streamlines[decomp_dim] = thisval
                        if 'df_streamlines' not in locals():
                            df_streamlines = tmp_streamlines.copy()
                        else:
                            df_streamlines = pd.concat([df_streamlines, tmp_streamlines], ignore_index=True)
                        # --------------------------------------------------

                        # reset n_samples
                        n_samples = n_samples_init  # probably unnecessary now that this is changed in a separate fxn, but can't hurt
                        # store wtd result
                        tau_bar_Ewtd_land[y,x,thisidx] = E_terrFrac * tau_wtd_mean
                        # moisture dist inland
                        moisture_dist_inland[y,x,thisidx] = dist_inland
                        # fraction of streamline over land 
                        land_frac_of_streamline[y,x,thisidx] = dist_inland / np.nanmax(dist)
                        # fraction of evap source coming from land ET 
                        land_frac_of_evapSource[y,x,thisidx] = land_Esource_frac

                        # repeat for decomposed ELs cases if necessary
                        if tau_ELs_decompose:
                            # [1] New E, old L and old flowfield
                            coast_step_E, nanindex_E, dist_E, dist_inland_E, tau_E, mu_E, wp_E, E0_E, P0_E, vq0_E, uq0_E, Fmag0_E, lfrac0_E, PminE_E, lat_save_E, lon_save_E = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                                                        Efit, Pfiti, vqfiti, uqfiti, Fmagfiti, lfracfit)
                            # compute E source fraction terrestrial
                            tau_wtd_mean_E, land_Esource_frac_E, E_terrFrac_E = fxn.terrestrial_E_frac(land_check, wp_E, coast_step_E, tau_E, E0_E, nanindex_E, n_samples)
                            # [2] New L, old E and old flowfield
                            coast_step_L, nanindex_L, dist_L, dist_inland_L, tau_L, mu_L, wp_L, E0_L, P0_L, vq0_L, uq0_L, Fmag0_L, lfrac0_L, PminE_L, lat_save_L, lon_save_L = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                                                        Efiti, Pfit, vqfiti, uqfiti, Fmagfit, lfracfit)
                            # compute E source fraction terrestrial
                            tau_wtd_mean_L, land_Esource_frac_L, E_terrFrac_L = fxn.terrestrial_E_frac(land_check, wp_L, coast_step_L, tau_L, E0_L, nanindex_L, n_samples)
                            # [3] New Path (s), old E and old L
                            coast_step_s, nanindex_s, dist_s, dist_inland_s, tau_s, mu_s, wp_s, E0_s, P0_s, vq0_s, uq0_s, Fmag0_s, lfrac0_s, PminE_s, lat_save_s, lon_save_s = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                                                        Efiti, Pfiti, vqfit, uqfit, Fmagfiti, lfracfit)
                            # compute E source fraction terrestrial
                            tau_wtd_mean_s, land_Esource_frac_s, E_terrFrac_s = fxn.terrestrial_E_frac(land_check, wp_s, coast_step_s, tau_s, E0_s, nanindex_s, n_samples)
                        
                                                
                        # resample to tau grid (instead of space) ----------
                        # (this is necessary to decompose E and L)
                        mu_taugrid, wp_taugrid, E0_taugrid, tausteps = fxn.dist_to_tau_coords(tau, nanindex, mu, wp, E0)

                        # TAU_BAR COMPUTATION
                        tau_bar[y,x,thisidx] = np.trapz(E0_taugrid * (1/mu_taugrid) * tausteps * math.e**(-tausteps)) / np.trapz(E0_taugrid * (1/mu_taugrid) * math.e**(-tausteps))
                        # get Dtau_bar
                        Dtau_bar[y,x,thisidx] = tau_bar[y,x,thisidx] - taubar_init

                        # tau decompositions -------------------------------
                        # decompose tau to local and regional components
                        # Note - unlike above, we integrate across space (not tau) so we can follow the spatial constraint on what's happening
                        if localEvapEff:
                            Dtau_bar_regionalEff_localEvap, Dtau_bar_localEff_localEvap = fxn.tau_spatial_decompose(local_threshold, x, y, thisidx,
                                                                                                dist, tau, E0, nanindex,
                                                                                                disti, taui, E0i, nanindexi,
                                                                                                Dtau_bar_regionalEff_localEvap, Dtau_bar_localEff_localEvap)
                        if Dtau_50_solve_localEvap:
                            Dtau_bar_50_dist_localEvap = fxn.dist_to_Dtaubar_frac_LocalEvap(x, y, thisidx, Dtau_bar_50_dist_localEvap, Dtau_bar, 
                                                                        dist, tau, E0, nanindex,
                                                                        disti, taui, E0i, nanindexi,
                                                                        Dtau_frac=Dtau_frac_upwindSolve, Dtau_bar_dist_percError = 1)
                        # get local v regional effects of changing mu and recompute tau upwind
                        if UpwindEff:  
                            tau_localConst, E0_localConst, tau_regionalConst, E0_regionalConst, Fmag0_localConst, Fmag0_regionalConst, P0_localConst, P0_regionalConst, dist_localConst, dist_regionalConst = fxn.local_plus_regional_tauStreamline(local_threshold, dx,
                                                                                                                                                                                                                                            P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
                                                                                                                                                                                                                                            P0, Fmag0, E0, dist, nanindex, lfrac0)
                            # get tau_bar due to regional change
                            taubar_localConst = np.trapz(tau_localConst*E0_localConst*math.e**(-tau_localConst))/np.trapz(E0_localConst*math.e**(-tau_localConst))
                            # and due to local change
                            taubar_regionalConst = np.trapz(tau_regionalConst*E0_regionalConst*math.e**(-tau_regionalConst))/np.trapz(E0_regionalConst*math.e**(-tau_regionalConst))
                            # update arrays
                            Dtau_bar_regionalEff_upwind[y,x,thisidx] = taubar_localConst - taubar_init 
                            Dtau_bar_localEff_upwind[y,x,thisidx] = taubar_regionalConst - taubar_init
                            # add streamline 
                            if save_upwindEff_streamline:
                                # construct the output dataframe ------------------
                                # set "df_streamlines" to none each time so it just gives the current streamline and we add decomp_dim and concatenate after the fact
                                no_streamlines_yet = True   # hold always true for this version
                                PminE_regionalConst = P0_regionalConst - E0_regionalConst
                                PminE_localConst = P0_localConst - E0_localConst
                                # additional terms
                                mu_regionalConst = P0_regionalConst / Fmag0_regionalConst
                                wp_regionalConst = (E0_regionalConst * math.e**(-tau_regionalConst)) / (np.trapz(E0_regionalConst * math.e**(-tau_regionalConst)))
                                mu_localConst = P0_localConst / Fmag0_localConst
                                wp_localConst = (E0_localConst * math.e**(-tau_localConst)) / (np.trapz(E0_localConst * math.e**(-tau_localConst)))
                                # define vars to ignore 
                                kwargs['streamline_ignore_vars'] = ["lat_save", "lon_save", "lfrac0", "nanindex"]
                                tmp_streamlines_localEff, _, _, _, _ = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                                    save_streamline, dist_regionalConst, lat_save, lon_save, E0_regionalConst, Fmag0_regionalConst, 
                                                                                                    P0_regionalConst, PminE_regionalConst, tau_regionalConst, mu_regionalConst, wp_regionalConst, lfrac0, nanindex, streamline_type='localChangeOnly',
                                                                                                    df_streamlines=None, **kwargs)
                                tmp_streamlines_regEff, _, _, _, _ = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                                    save_streamline, dist_localConst, lat_save, lon_save, E0_localConst, Fmag0_localConst, 
                                                                                                    P0_localConst, PminE_localConst, tau_localConst, mu_localConst, wp_localConst, lfrac0, nanindex, streamline_type='regionalChangeOnly',
                                                                                                    df_streamlines=None, **kwargs)
                                tmp_streamlines = pd.concat([tmp_streamlines_localEff, tmp_streamlines_regEff])
                                tmp_streamlines[decomp_dim] = thisval
                                if 'df_streamlines' not in locals():
                                    df_streamlines = tmp_streamlines.copy()
                                else:
                                    df_streamlines = pd.concat([df_streamlines, tmp_streamlines], ignore_index=True)
                                kwargs['streamline_ignore_vars'] = None
                                # --------------------------------------------------

                        if Dtau_50_solve_strmline:
                            Dtau_bar_50_dist_strmline = fxn.dist_to_Dtaubar_frac_streamline(x, y, thisidx, dx, Dtau_bar, taubar_init,
                                                                                        P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
                                                                                        P0, Fmag0, E0, dist, nanindex, lfrac0,
                                                                                        Dtau_bar_50_dist_strmline,
                                                                                        Dtau_frac=Dtau_frac_upwindSolve, Dtau_bar_dist_percError = 1, troubleshoot=False)
                            
                        # decompose E, L, and path
                        if tau_ELs_decompose:
                            Dtau_bar_newE, Dtau_bar_newL, Dtau_bar_newPath = fxn.tau_decompose_ELPath(x, y, thisidx, taubar_init,
                                                                                                    E0_E, E0_L, E0_s,
                                                                                                    mu_E, mu_L, mu_s,
                                                                                                    wp_E, wp_L, wp_s,
                                                                                                    tau_E, tau_L, tau_s, 
                                                                                                    nanindex_E, nanindex_L, nanindex_s,
                                                                                                    Dtau_bar_newE, Dtau_bar_newL, Dtau_bar_newPath)
                
                        # track progress -------------------------------------------------------------
                        complete_percent = round((readout_counter / n_cells)*100,2)
                        l1 = "Run: " + this_runName  
                        l2 = "Step: " + str(this_step) + " of " + str(total_steps) 
                        l3 = "Step progress: " + str(complete_percent) + "% ... " 
                        l4 = "--------------------------------------------"
                        print(l1, l2, l3, l4, sep=os.linesep)
                        # ----------------------------------------------------------------------------
                        
                        # adjust counter
                        readout_counter += 1

                    # end decomp_dim loop
                    idx_counter += 1

            # bring arrays together in dict object
            all_arrays = {"tau_bar": tau_bar, "tau_bar_Ewtd_land": tau_bar_Ewtd_land,
                                "moisture_dist_inland": moisture_dist_inland, 
                                "land_frac_of_streamline": land_frac_of_streamline,
                                "land_frac_of_evapSource": land_frac_of_evapSource,
                                }
            if localEvapEff:
                all_arrays["Dtau_bar_regionalEff_localEvap"] = Dtau_bar_regionalEff_localEvap
                all_arrays["Dtau_bar_localEff_localEvap"] = Dtau_bar_localEff_localEvap
            if UpwindEff:
                all_arrays["Dtau_bar_regionalEff_upwind"] = Dtau_bar_regionalEff_upwind
                all_arrays["Dtau_bar_localEff_upwind"] = Dtau_bar_localEff_upwind
            if tau_ELs_decompose:
                all_arrays["Dtau_bar_newE"] = Dtau_bar_newE
                all_arrays["Dtau_bar_newL"] = Dtau_bar_newL
                all_arrays["Dtau_bar_newPath"] = Dtau_bar_newPath
            if localEvapEff or UpwindEff or tau_ELs_decompose:
                all_arrays["Dtau_bar"] = Dtau_bar
            if Dtau_50_solve_localEvap:
                all_arrays["Dtau_bar_50_dist_localEvap"] = Dtau_bar_50_dist_localEvap
            if Dtau_50_solve_strmline:
                all_arrays["Dtau_bar_50_dist_upwind"] = Dtau_bar_50_dist_strmline

            # restore to original resolution
            if res_adj:
                old_Y, old_X = np.meshgrid(lon, lat)
                new_y = ds.lat.values
                new_x = ds.lon.values
                new_Y, new_X = np.meshgrid(new_x, new_y)
                # sub arrays to tmp dict
                coarse_arrays = all_arrays.copy()  # store
                all_arrays = {}  # make empty to refill in the loop
                # loop through all arrays
                for name, arr in coarse_arrays.items():
                    # initialize new array
                    new_arr = np.full((ds.lat.shape[0],ds.lon.shape[0], len(kwargs.get('decomp'))), np.nan)

                    # loop through decomp_dim
                    for decompdim in range(len(ds[decomp_dim])):
                        tmp_arr = arr[:,:,decompdim]
                        new_arr[:,:,decompdim] = griddata((old_Y.flatten(), old_X.flatten()), tmp_arr.flatten(), (new_Y, new_X), method='linear')
                    # store in all_arr
                    all_arrays[name] = new_arr

        
        
        # COORD-LIST LOOP -----------------------------------------------------------------
        if compute_method == 'coord_list':   
            decomp_slice = range(len(ds[decomp_dim])-1)  # subtract 1 because we compute base case outside the decomposition dim loop
            case_idxs = np.where(ds[decomp_dim].values != baseline)[0]   # all indices except baseline  
            idx_counter = 1
            readout_counter = 1
            n_cells = len(iso_coords)
            n_iters = n_cells * len(decomp_slice)
            # -- LOOP THROUGH COORDS -- 
            for coord in range(n_cells):
                # coordinates at point of precipitation
                lat0 = iso_coords['lat'].iloc[coord]
                lon0 = iso_coords['lon'].iloc[coord]

                # check whether to save streamline
                if streamline_collect == "some":
                    if iso_coords['collect_streamline'].iloc[coord] == "Y":
                        save_streamline = True
                elif streamline_collect == "all":
                    save_streamline = True

                # convert to 0-360 if necessary
                if lon0 < 0:
                    lon0 = lon0 % 360

                # get x and y indices of nearest lat/lon in grid (for storing result)
                latval = ds.sel(lat=lat0, method='nearest').lat.values
                lonval = ds.sel(lon=lon0, method='nearest').lon.values
                x = np.where(ds.lon.values == lonval)[0][0]
                y = np.where(ds.lat.values == latval)[0][0]

                # check if it's on land
                land_check = ds['LANDFRAC'].sel(lat=lat0, lon=lon0, method='nearest').values[0]
                if land_only and land_check == 0.:
                    moisture_dist_inland[y,x,:] = 0
                    # track progress -------------------------------------------------------------
                    complete_percent = round((idx_counter / n_iters)*100,2)
                    l1 = "...skipping tau over the ocean " + " -- " + str(complete_percent) + "%" + " complete"
                    l2 = "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ "
                    print(l1, l2, sep=os.linesep)
                    # ----------------------------------------------------------------------------

                    idx_counter += 1 # move to next index
                    continue
                elif land_check < 0.6:
                    moisture_dist_inland[y,x,baseline_idx] = 0

                # compute tau and climatology across streamline
                coast_stepi, nanindexi, disti, dist_inlandi, taui, mui, wpi, E0i, P0i, vq0i, uq0i, Fmag0i, lfrac0i, PminEi, lat_savei, lon_savei = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                            Efiti, Pfiti, vqfiti, uqfiti, Fmagfiti, lfracfiti)
                # compute E source fraction terrestrial
                tau_wtd_meani, land_Esource_fraci, E_terrFraci = fxn.terrestrial_E_frac(land_check, wpi, coast_stepi, taui, E0i, nanindexi, n_samples)
                
                # construct the output dataframe ------------------
                # set "df_streamlines" to none each time so it just gives the current streamline and we add decomp_dim and concatenate after the fact
                no_streamlines_yet = True   # hold always true for this version
                df_streamlinesi, _, _, _, _ = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                    save_streamline, disti, lat_savei, lon_savei, E0i, Fmag0i, P0i, PminEi, 
                                                                                    taui, mui, wpi, lfrac0i, nanindexi, streamline_type='full_model',
                                                                                    df_streamlines=None, **kwargs)
                df_streamlinesi[decomp_dim] = baseline
                if 'df_streamlines' not in locals():
                    df_streamlines = df_streamlinesi.copy()
                else:
                    df_streamlines = pd.concat([df_streamlines, df_streamlinesi], ignore_index=True)
                # --------------------------------------------------

                # reset n_samples
                n_samples = n_samples_init  # probably unnecessary now that this is changed in a separate fxn, but can't hurt
                # store wtd result
                tau_bar_Ewtd_land[y,x,baseline_idx] = E_terrFraci * tau_wtd_meani
                # moisture dist inland
                moisture_dist_inland[y,x,baseline_idx] = dist_inlandi
                # fraction of streamline over land 
                land_frac_of_streamline[y,x,baseline_idx] = dist_inlandi / np.nanmax(disti)
                # fraction of evap source coming from land ET 
                land_frac_of_evapSource[y,x,baseline_idx] = land_Esource_fraci

                # resample to tau grid (instead of space) ----------
                # (this is necessary to decompose E and L)
                mu_taugridi, wp_taugridi, E0_taugridi, taustepsi = fxn.dist_to_tau_coords(taui, nanindexi, mui, wpi, E0i)

                # SOLVE TAU-BAR
                tau_bar[y,x,baseline_idx] = np.trapz(E0_taugridi * (1/mu_taugridi) * taustepsi * math.e**(-taustepsi)) / np.trapz(E0_taugridi * (1/mu_taugridi) * math.e**(-taustepsi))
                taubar_init = tau_bar[y,x,baseline_idx] 

                # -- LOOP THROUGH DECOMPOSITION DIM -- 
                for ddim in decomp_slice:
                    thisidx = case_idxs[ddim]
                    thisval = ds[decomp_dim].values[thisidx]
                    # build the grid interpolator 
                    # get interpolation functions for baseline case 
                    Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit = fxn.tau_interpolatorInitialize(ds.lon.values, ds.lat.values, E[:,:,thisidx], P[:,:,thisidx], 
                                                                                                UQ[:,:,thisidx], VQ[:,:,thisidx], lfrac[:,:,thisidx])
                    # compute tau and climatology across streamline
                    coast_step, nanindex, dist, dist_inland, tau, mu, wp, E0, P0, vq0, uq0, Fmag0, lfrac0, PminE, lat_save, lon_save = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                Efit, Pfit, vqfit, uqfit, Fmagfit, lfracfit)
                    # compute E source fraction terrestrial
                    tau_wtd_mean, land_Esource_frac, E_terrFrac = fxn.terrestrial_E_frac(land_check, wp, coast_step, tau, E0, nanindex, n_samples)
                    
                    # update streamlines
                    # construct the output dataframe ------------------
                    # set "df_streamlines" to none each time so it just gives the current streamline and we add decomp_dim and concatenate after the fact
                    no_streamlines_yet = True   # hold always true for this version
                    tmp_streamlines, _, _, _, _ = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                        save_streamline, dist, lat_save, lon_save, E0, Fmag0, P0, PminE, 
                                                                                        tau, mu, wp, lfrac0, nanindex, streamline_type='full_model',
                                                                                        df_streamlines=None, **kwargs)
                    tmp_streamlines[decomp_dim] = thisval
                    if 'df_streamlines' not in locals():
                        df_streamlines = tmp_streamlines.copy()
                    else:
                        df_streamlines = pd.concat([df_streamlines, tmp_streamlines], ignore_index=True)
                    # --------------------------------------------------

                    # reset n_samples
                    n_samples = n_samples_init  # probably unnecessary now that this is changed in a separate fxn, but can't hurt
                    # store wtd result
                    tau_bar_Ewtd_land[y,x,thisidx] = E_terrFrac * tau_wtd_mean
                    # moisture dist inland
                    moisture_dist_inland[y,x,thisidx] = dist_inland
                    # fraction of streamline over land 
                    land_frac_of_streamline[y,x,thisidx] = dist_inland / np.nanmax(dist)
                    # fraction of evap source coming from land ET 
                    land_frac_of_evapSource[y,x,thisidx] = land_Esource_frac

                    # repeat for decomposed ELs cases if necessary
                    if tau_ELs_decompose:
                        # [1] New E, old L and old flowfield
                        coast_step_E, nanindex_E, dist_E, dist_inland_E, tau_E, mu_E, wp_E, E0_E, P0_E, vq0_E, uq0_E, Fmag0_E, lfrac0_E, PminE_E, lat_save_E, lon_save_E = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                                                    Efit, Pfiti, vqfiti, uqfiti, Fmagfiti, lfracfit)
                        # compute E source fraction terrestrial
                        tau_wtd_mean_E, land_Esource_frac_E, E_terrFrac_E = fxn.terrestrial_E_frac(land_check, wp_E, coast_step_E, tau_E, E0_E, nanindex_E, n_samples)
                        # [2] New L, old E and old flowfield
                        coast_step_L, nanindex_L, dist_L, dist_inland_L, tau_L, mu_L, wp_L, E0_L, P0_L, vq0_L, uq0_L, Fmag0_L, lfrac0_L, PminE_L, lat_save_L, lon_save_L = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                                                    Efiti, Pfit, vqfiti, uqfiti, Fmagfit, lfracfit)
                        # compute E source fraction terrestrial
                        tau_wtd_mean_L, land_Esource_frac_L, E_terrFrac_L = fxn.terrestrial_E_frac(land_check, wp_L, coast_step_L, tau_L, E0_L, nanindex_L, n_samples)
                        # [3] New Path (s), old E and old L
                        coast_step_s, nanindex_s, dist_s, dist_inland_s, tau_s, mu_s, wp_s, E0_s, P0_s, vq0_s, uq0_s, Fmag0_s, lfrac0_s, PminE_s, lat_save_s, lon_save_s = fxn.tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, 
                                                                                                                                                                                                    Efiti, Pfiti, vqfit, uqfit, Fmagfiti, lfracfit)
                        # compute E source fraction terrestrial
                        tau_wtd_mean_s, land_Esource_frac_s, E_terrFrac_s = fxn.terrestrial_E_frac(land_check, wp_s, coast_step_s, tau_s, E0_s, nanindex_s, n_samples)

                    # resample to tau grid (instead of space) ----------
                    # (this is necessary to decompose E and L)
                    mu_taugrid, wp_taugrid, E0_taugrid, tausteps = fxn.dist_to_tau_coords(tau, nanindex, mu, wp, E0)

                    # TAU_BAR COMPUTATION
                    tau_bar[y,x,thisidx] = np.trapz(E0_taugrid * (1/mu_taugrid) * tausteps * math.e**(-tausteps)) / np.trapz(E0_taugrid * (1/mu_taugrid) * math.e**(-tausteps))
                    # get Dtau_bar
                    Dtau_bar[y,x,thisidx] = tau_bar[y,x,thisidx] - taubar_init

                     # tau decompositions -------------------------------
                    # decompose tau to local and regional components
                    # Note - unlike above, we integrate across space (not tau) so we can follow the spatial constraint on what's happening
                    if localEvapEff:
                        Dtau_bar_regionalEff_localEvap, Dtau_bar_localEff_localEvap = fxn.tau_spatial_decompose(local_threshold, x, y, thisidx,
                                                                                            dist, tau, E0, nanindex,
                                                                                            disti, taui, E0i, nanindexi,
                                                                                            Dtau_bar_regionalEff_localEvap, Dtau_bar_localEff_localEvap)
                    if Dtau_50_solve_localEvap:
                        Dtau_bar_50_dist_localEvap = fxn.dist_to_Dtaubar_frac_LocalEvap(x, y, thisidx, Dtau_bar_50_dist_localEvap, Dtau_bar, 
                                                                    dist, tau, E0, nanindex,
                                                                    disti, taui, E0i, nanindexi,
                                                                    Dtau_frac=Dtau_frac_upwindSolve, Dtau_bar_dist_percError = 1)
                    # get local v regional effects of changing mu and recompute tau upwind
                    if UpwindEff:  
                        tau_localConst, E0_localConst, tau_regionalConst, E0_regionalConst, Fmag0_localConst, Fmag0_regionalConst, P0_localConst, P0_regionalConst, dist_localConst, dist_regionalConst = fxn.local_plus_regional_tauStreamline(local_threshold, dx,
                                                                                                                                                                                                                                                    P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
                                                                                                                                                                                                                                                    P0, Fmag0, E0, dist, nanindex, lfrac0)
                        # get tau_bar due to regional change
                        taubar_localConst = np.trapz(tau_localConst*E0_localConst*math.e**(-tau_localConst))/np.trapz(E0_localConst*math.e**(-tau_localConst))
                        # and due to local change
                        taubar_regionalConst = np.trapz(tau_regionalConst*E0_regionalConst*math.e**(-tau_regionalConst))/np.trapz(E0_regionalConst*math.e**(-tau_regionalConst))
                        # update arrays
                        Dtau_bar_regionalEff_upwind[y,x,thisidx] = taubar_localConst - taubar_init 
                        Dtau_bar_localEff_upwind[y,x,thisidx] = taubar_regionalConst - taubar_init
                        # add streamline 
                        if save_upwindEff_streamline:
                            # construct the output dataframe ------------------
                            # set "df_streamlines" to none each time so it just gives the current streamline and we add decomp_dim and concatenate after the fact
                            no_streamlines_yet = True   # hold always true for this version
                            PminE_regionalConst = P0_regionalConst - E0_regionalConst
                            PminE_localConst = P0_localConst - E0_localConst
                            # additional terms
                            mu_regionalConst = P0_regionalConst / Fmag0_regionalConst
                            wp_regionalConst = (E0_regionalConst * math.e**(-tau_regionalConst)) / (np.trapz(E0_regionalConst * math.e**(-tau_regionalConst)))
                            mu_localConst = P0_localConst / Fmag0_localConst
                            wp_localConst = (E0_localConst * math.e**(-tau_localConst)) / (np.trapz(E0_localConst * math.e**(-tau_localConst)))
                            # define vars to ignore 
                            kwargs['streamline_ignore_vars'] = ["lat_save", "lon_save", "lfrac0", "nanindex"]
                            tmp_streamlines_localEff, _, _, _, _ = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                                save_streamline, dist_regionalConst, lat_save, lon_save, E0_regionalConst, Fmag0_regionalConst, 
                                                                                                P0_regionalConst, PminE_regionalConst, tau_regionalConst, mu_regionalConst, wp_regionalConst, lfrac0, nanindex, streamline_type='localChangeOnly',
                                                                                                df_streamlines=None, **kwargs)
                            tmp_streamlines_regEff, _, _, _, _ = fxn.build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                                                                                                save_streamline, dist_localConst, lat_save, lon_save, E0_localConst, Fmag0_localConst, 
                                                                                                P0_localConst, PminE_localConst, tau_localConst, mu_localConst, wp_localConst, lfrac0, nanindex, streamline_type='regionalChangeOnly',
                                                                                                df_streamlines=None, **kwargs)
                            tmp_streamlines = pd.concat([tmp_streamlines_localEff, tmp_streamlines_regEff])
                            tmp_streamlines[decomp_dim] = thisval
                            if 'df_streamlines' not in locals():
                                df_streamlines = tmp_streamlines.copy()
                            else:
                                df_streamlines = pd.concat([df_streamlines, tmp_streamlines], ignore_index=True)
                            # reset ignore vars to none
                            kwargs['streamline_ignore_vars'] = None
                            # --------------------------------------------------


                    if Dtau_50_solve_strmline:
                        Dtau_bar_50_dist_strmline = fxn.dist_to_Dtaubar_frac_streamline(x, y, thisidx, dx, Dtau_bar, taubar_init,
                                                                                    P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
                                                                                    P0, Fmag0, E0, dist, nanindex, lfrac0,
                                                                                    Dtau_bar_50_dist_strmline,
                                                                                    Dtau_frac=Dtau_frac_upwindSolve, Dtau_bar_dist_percError = 1, troubleshoot=False)
                        
                    # decompose E, L, and path
                    if tau_ELs_decompose:
                        Dtau_bar_newE, Dtau_bar_newL, Dtau_bar_newPath = fxn.tau_decompose_ELPath(x, y, thisidx, taubar_init,
                                                                                                E0_E, E0_L, E0_s,
                                                                                                mu_E, mu_L, mu_s,
                                                                                                wp_E, wp_L, wp_s,
                                                                                                tau_E, tau_L, tau_s, 
                                                                                                nanindex_E, nanindex_L, nanindex_s,
                                                                                                Dtau_bar_newE, Dtau_bar_newL, Dtau_bar_newPath)
            
                    # track progress -------------------------------------------------------------
                    complete_percent = round((readout_counter / n_iters)*100,2)
                    l1 = "Year: " + year_printout + ", Slice: " + slice_printout + " -- " + str(complete_percent) + "%" + " complete"
                    l2 = "(" + yrslc_remain_printout + " Year-Slices remaining)"
                    l3 = "--------------------------------------------"
                    print(l1, l2, l3, sep=os.linesep)
                    # ----------------------------------------------------------------------------
                    
                    # adjust counter
                    readout_counter += 1

                # end decomp_dim loop
                idx_counter += 1

            # bring arrays together in dict object
            all_arrays = {"tau_bar": tau_bar, "tau_bar_Ewtd_land": tau_bar_Ewtd_land,
                                "moisture_dist_inland": moisture_dist_inland, 
                                "land_frac_of_streamline": land_frac_of_streamline,
                                "land_frac_of_evapSource": land_frac_of_evapSource,
                                }
            if localEvapEff:
                all_arrays["Dtau_bar_regionalEff_localEvap"] = Dtau_bar_regionalEff_localEvap
                all_arrays["Dtau_bar_localEff_localEvap"] = Dtau_bar_localEff_localEvap
            if UpwindEff:
                all_arrays["Dtau_bar_regionalEff_upwind"] = Dtau_bar_regionalEff_upwind
                all_arrays["Dtau_bar_localEff_upwind"] = Dtau_bar_localEff_upwind
            if tau_ELs_decompose:
                all_arrays["Dtau_bar_newE"] = Dtau_bar_newE
                all_arrays["Dtau_bar_newL"] = Dtau_bar_newL
                all_arrays["Dtau_bar_newPath"] = Dtau_bar_newPath
            if localEvapEff or UpwindEff or tau_ELs_decompose:
                all_arrays["Dtau_bar"] = Dtau_bar
            if Dtau_50_solve_localEvap:
                all_arrays["Dtau_bar_50_dist_localEvap"] = Dtau_bar_50_dist_localEvap
            if Dtau_50_solve_strmline:
                all_arrays["Dtau_bar_50_dist_upwind"] = Dtau_bar_50_dist_strmline


        # ADD RESULT TO DS
        tau_name = "tau_bar"
        moisture_dist_name = "moisture_dist_inland"
        tau_bar_EWL_name = "tau_bar_wtdLandEvap"
        landfrac_streamline_name = "streamline_frac_land"
        land_frac_of_evapSource_name = "Esource_frac_land"
        # get to data array
        ds[tau_name] = xr.DataArray(all_arrays["tau_bar"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=tau_name)
        ds[moisture_dist_name] = xr.DataArray(all_arrays["moisture_dist_inland"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=moisture_dist_name)
        ds[tau_bar_EWL_name] = xr.DataArray(all_arrays["tau_bar_Ewtd_land"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=tau_bar_EWL_name)
        ds[landfrac_streamline_name] = xr.DataArray(all_arrays["land_frac_of_streamline"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=landfrac_streamline_name)
        ds[land_frac_of_evapSource_name] = xr.DataArray(all_arrays["land_frac_of_evapSource"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=land_frac_of_evapSource_name)
        
        if localEvapEff:
            # get sum 
            Dtau_bar_localRegional_localEvap_sum = all_arrays["Dtau_bar_localEff_localEvap"] + all_arrays["Dtau_bar_regionalEff_localEvap"]
            # assign names
            Dtau_leff_name = "Dtau_bar_localEffect_localEvap_l" + str(np.around(local_threshold))
            Dtau_reff_name = "Dtau_bar_regEffect_localEvap_g" + str(np.around(local_threshold))
            Dtau_lrsum_name = "Dtau_bar_localRegional_localEvap_effectSum" 
            # add to ds
            ds[Dtau_leff_name] = xr.DataArray(all_arrays["Dtau_bar_localEff_localEvap"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_leff_name)
            ds[Dtau_reff_name] = xr.DataArray(all_arrays["Dtau_bar_regionalEff_localEvap"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_reff_name)
            ds[Dtau_lrsum_name] = xr.DataArray(Dtau_bar_localRegional_localEvap_sum[:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_lrsum_name)
        if UpwindEff:
            # get sum 
            Dtau_bar_localRegional_upwind_sum = all_arrays["Dtau_bar_localEff_upwind"] + all_arrays["Dtau_bar_regionalEff_upwind"]
            # assign names
            Dtau_leff_upwind_name = "Dtau_bar_localEffect_upwind_l" + str(np.around(local_threshold))
            Dtau_reff_upwind_name = "Dtau_bar_regEffect_upwind_g" + str(np.around(local_threshold))
            Dtau_lrsum_upwind_name = "Dtau_bar_localRegional_upwind_effectSum" 
            # add to ds
            ds[Dtau_leff_upwind_name] = xr.DataArray(all_arrays["Dtau_bar_localEff_upwind"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_leff_upwind_name)
            ds[Dtau_reff_upwind_name] = xr.DataArray(all_arrays["Dtau_bar_regionalEff_upwind"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_reff_upwind_name)
            ds[Dtau_lrsum_upwind_name] = xr.DataArray(Dtau_bar_localRegional_upwind_sum[:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_lrsum_upwind_name)
        if tau_ELs_decompose:
            # get sum 
            Dtau_bar_ELs_sum = all_arrays["Dtau_bar_newE"] + all_arrays["Dtau_bar_newL"] + all_arrays["Dtau_bar_newPath"]
            # assign names
            Dtau_Ename = "Dtau_bar_Eeffect"
            Dtau_Lname = "Dtau_bar_Leffect"
            Dtau_sname = "Dtau_bar_PathEffect"
            Dtau_ELs_sumname = "Dtau_bar_ELs_effectSum"
            # add to ds
            ds[Dtau_Ename] = xr.DataArray(all_arrays["Dtau_bar_newE"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_Ename)
            ds[Dtau_Lname] = xr.DataArray(all_arrays["Dtau_bar_newL"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_Lname)
            ds[Dtau_sname] = xr.DataArray(all_arrays["Dtau_bar_newPath"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_sname)
            ds[Dtau_ELs_sumname] = xr.DataArray(Dtau_bar_ELs_sum[:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_ELs_sumname)
        if Dtau_50_solve_localEvap:
            # name 
            Dtau_50_name_localEvap = "Dtau_bar_distTo"+str(int(np.around(Dtau_frac_upwindSolve*100,0)))+"perc_localEvap"
            ds[Dtau_50_name_localEvap] = xr.DataArray(all_arrays["Dtau_bar_50_dist_localEvap"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_50_name_localEvap)
        if Dtau_50_solve_strmline:
            # name 
            Dtau_50_name_upwind = "Dtau_bar_distTo"+str(int(np.around(Dtau_frac_upwindSolve*100,0)))+"perc_upwind"
            ds[Dtau_50_name_upwind] = xr.DataArray(all_arrays["Dtau_bar_50_dist_upwind"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_50_name_upwind)        
        if tau_ELs_decompose or localEvapEff or UpwindEff:
            Dtau_name = "Dtau_bar"
            ds[Dtau_name] = xr.DataArray(all_arrays["Dtau_bar"][:,:,:], 
                                    dims=('lat', 'lon', decomp_dim),
                                    coords={'lat': ds.lat.values, 'lon': ds.lon.values, decomp_dim: ds[decomp_dim].values},
                                    name=Dtau_name)

        # return result
        return ds, df_streamlines
        




            

# ----------------------------------------------------------------
class Run:

    # [1] READ IN INPUTS AND INITIALIZE
    @staticmethod
    def inputfiles(**kwargs):
        # read in clim and forcing ------------------------------
        in_path = os.path.join(kwargs.get('run_path'), kwargs.get('run_name'), kwargs.get('input_dir'))
        in_clim_path = os.path.join(in_path, kwargs.get('clim_fn'))
        # check if force name is none
        if kwargs.get('force_fn') is None:
            kwargs['force_fn'] = 'None'  # (otherwise we get an error when checking path... lazy but simple solution :)
        in_force_path = os.path.join(in_path, kwargs.get('force_fn'))
        # check if clim file exist
        if os.path.exists(in_clim_path):
            # bring in clim
            ds_clim_all = xr.open_dataset(in_clim_path)
        else:
            raise ValueError("I can't find a climate input file, check that the filename is correct in initialize.py and that it's located in the /inputs directory")
        # check for forcing file 
        if os.path.exists(in_force_path):
            ds_force_all = xr.open_dataset(in_force_path)
        else:
            ds_force_all = None
            warnings.warn("Couldn't find a forcing file, I'll run just the isotope module (no EBM)", UserWarning)
        
        # read in coords if necessary ----------------------------
        # read in coordinates if necessary
        if kwargs.get('Format_for_where_to_compute_tau-bar') == "coord_list":
            coord_fn = os.path.join(in_path, kwargs['CoordList_filename'])
            df_coords = pd.read_csv(coord_fn)
            # remove any rows w/ missing data
            notmissing = df_coords[['lat', 'lon']].notnull().all(axis=1)
            df_coords = df_coords[notmissing]

            if kwargs.get('Collect_StreamlineDat_For-') == "some":
                if "collect_streamline" not in df_coords.columns:
                    raise ValueError("If you only want me to output 'some' isotope streamline data, you need a 'collect_streamline' column in isotope coord_list.csv")

        if kwargs.get('Format_for_where_to_compute_tau-bar') == "bbox":
            if kwargs.get('Collect_StreamlineDat_For-') == "some":
                coord_fn = os.path.join(in_path, kwargs['CoordList_filename'])
                df_coords = pd.read_csv(coord_fn)

                if "collect_streamline" not in df_coords.columns:
                    raise ValueError("If you only want me to output 'some' isotope streamline data, you need a 'collect_streamline' column in isotope coord_list.csv")

        if 'df_coords' not in locals():
            df_coords = None

        # determine runtype --------------------------------------
        if ds_force_all is not None:
            if kwargs.get('SolveIsotopes'):
                runtype = "full_model"
            else:
                runtype = "EBM_noIsotopes"
        else:
            runtype = "isotopes_noEBM"
        # add back to kwargs
        kwargs['runtype'] = runtype

        # return results
        return ds_clim_all, ds_force_all, df_coords, kwargs


    # [2] SOLVE JUST CLIMATOLOGY 
    @staticmethod
    def ebmsolve(ds_clim_all, ds_force_all, slice_steps, tsteps, **kwargs):
        # here, no matter the runtype, we'll solve climatology (which EBM functions are 
        # used is filtered by this function itself)
        rtype_init = kwargs.get('runtype')   # we'll add it back in later
        kwargs['runtype'] = "full_model"  # as long as it doesn't say "noEBM", you're good

        # if no ds_force, then only solve ebm -----------------------------
        if ds_force_all is None:  # (no forcing, just solve energetics for each case)
            # do we need to loop through time?
            if len(tsteps) <= 1: # -- no time loop
                for slc in slice_steps:
                    ds_clim = ds_clim_all.sel(yr_slice = slc)
                    # hydroclim partitioning (not necessary, but why not)
                    # hydroclim partitioning 
                    ds_clim = Hydroclim.orog_partition(ds_clim, apply_to_var="PRECT", **kwargs)
                    ds_clim = Hydroclim.orog_partition(ds_clim, apply_to_var="Q", **kwargs)
                    # get energetics
                    # ds_clim = EBM.divMSE(ds_clim)
                    # efe
                    # ds_clim = EBM.EFE_calc_heaviside(ds_clim, **kwargs)
                    # and efpm
                    # contours_filtered1 = EBM.stitch_efpm_contours(ds_clim, **kwargs)
                    # run this out of the if tstep==0 b/c we overwrite df1 later (so we need to re-overwrite it here)
                    # df1 = EBM.efpm_rank(contours=contours_filtered1, dataset=ds_clim, case_or_ctrl='ctrl', **kwargs)
                    # ds_clim = EBM.efpm_contours_to_mask(ds_clim, df1, **kwargs)

                    # bring together
                    if slc == slice_steps[0]: # combine control and first case scenario
                        outds = ds_clim.assign_coords(time = 0.)
                        outds = outds.expand_dims(dim='time')
                        outds = outds.expand_dims(yr_slice = [slc], axis=3)
                    # merge result
                    else:
                        tds = ds_clim.assign_coords(time = 0.)
                        tds = tds.expand_dims(dim='time')
                        tds = tds.expand_dims(yr_slice = [slc], axis=3)
                        outds = xr.merge([outds, tds])
                        
            else: # there is a time loop here... 
                for slc in slice_steps:
                    for tstep in tsteps:  
                        ds_clim = ds_clim_all.sel(yr_slice = slc, time=tstep)
                        # hydroclim partitioning (not necessary, but why not)
                        # hydroclim partitioning 
                        ds_clim = Hydroclim.orog_partition(ds_clim, apply_to_var="PRECT", **kwargs)
                        ds_clim = Hydroclim.orog_partition(ds_clim, apply_to_var="Q", **kwargs)
                        # get energetics
                        # ds_clim = EBM.divMSE(ds_clim)
                        # efe
                        # ds_clim = EBM.EFE_calc_heaviside(ds_clim, **kwargs)
                        # and efpm
                        # contours_filtered1 = EBM.stitch_efpm_contours(ds_clim, **kwargs)
                        # run this out of the if tstep==0 b/c we overwrite df1 later (so we need to re-overwrite it here)
                        # df1 = EBM.efpm_rank(contours=contours_filtered1, dataset=ds_clim, case_or_ctrl='ctrl', **kwargs)
                        # ds_clim = EBM.efpm_contours_to_mask(ds_clim, df1, **kwargs)

                        # bring together
                        if tstep == 0 and slc == slice_steps[0]: # combine control and first case scenario
                            outds = ds_clim.assign_coords(time=ds_clim_all.time[tstep])
                            outds = outds.expand_dims(dim='time')
                            outds = outds.expand_dims(yr_slice = [slc], axis=3)
                        # merge result
                        else:
                            tds = ds_clim.assign_coords(time=ds_clim_all.time[tstep])
                            tds = tds.expand_dims(dim='time')
                            tds = tds.expand_dims(yr_slice = [slc], axis=3)
                            outds = xr.merge([outds, tds])
        
        else: # here we have a forcing file and we assume it has a time component
                # ... SAVE TIME STEP RESULT ----------- 
                # assign time coord back
                if 'outds' not in locals():
                    outds = outds.expand_dims(yr_slice = [slc], axis=3)
                else:
                    tds = tds.expand_dims(yr_slice = [slc], axis=3)
                    # and merge
                    outds = xr.merge([outds, tds])                  

           
        # replace runtype 
        kwargs['runtype'] = rtype_init
        
        # return result
        return outds



    # [3] LOOP THROUGH 
    @staticmethod
    def loop_ebm_iso(ds_clim_all, ds_force_all, slice_steps, tsteps, 
                    spatial_decomposition, clim_decomposition, tau_decompByTstep, 
                    df_coords=None, **kwargs):
        # do we need to decompose tau at all?
        if spatial_decomposition or clim_decomposition:
            tau_decompose = True
            if tau_decompByTstep:
                decompose_dim, loop_dim = "time", "yr_slice"
            else:
                decompose_dim, loop_dim = "yr_slice", "time"
        else: 
            tau_decompose = False
            decompose_dim = None
        # ------------------------------------
        
        # first solve climatology ------------
        ds_ebm = Run.ebmsolve(ds_clim_all, ds_force_all, slice_steps, tsteps, **kwargs)

        # now run through isotopes -----------
        if not tau_decompose:   # with NO tau decomposition
            # update tstep based on EBM (orig. tstep is just n forcing steps, doesn't include time=1)
            ebm_tsteps = range(len(ds_ebm.time.values))
            kwargs['yr_slices_remaining_int'] = len(slice_steps) * len(ds_ebm.time.values) + 1  # add one b/c we subtract at start of tstep to keep all trackers together
            kwargs['yr_slices_remaining'] = str(kwargs['yr_slices_remaining_int'])
            # loop over time
            for tstep in ebm_tsteps:
                if 'time' in ds_ebm.dims:
                    ds_tstep = ds_ebm.isel(time=tstep)
                    this_timestep = ds_ebm.time.values[tstep]
                    kwargs['ts_year'] = str(int(this_timestep))
                else:
                    ds_tstep = ds_ebm.copy()
                    this_timestep = 0
                    kwargs['ts_year'] = str(int(this_timestep))
                
                # loop over slice
                for slc in slice_steps:
                    # ... update timing
                    kwargs['slice_of_year'] = slc
                    kwargs['yr_slices_remaining_int'] -= 1
                    kwargs['yr_slices_remaining'] = str(kwargs['yr_slices_remaining_int'])
                    # --------------------------------
                    ds_in = ds_tstep.sel(yr_slice=slc)
                    # solve for tau
                    ds_in, df_streamlines = Isotopes.taubar_noDecomposition(ds_in, iso_coords=df_coords, **kwargs)
                    # add info to df_streamlines
                    df_streamlines['yr_slice'] = slc
                    df_streamlines['time'] = this_timestep

                    # create output data
                    # ... SAVE SLICE RESULT ----------- 
                    ds_out = ds_in.assign_coords(time=this_timestep)
                    ds_out = ds_out.expand_dims(dim='time')
                    ds_out = ds_out.expand_dims(yr_slice = [slc], axis=3)
                    if tstep == 0 and slc == slice_steps[0]:
                        outds = ds_out.copy()
                        outstreamlines = df_streamlines
                    else:
                        outds = xr.merge([outds, ds_out])
                        outstreamlines = pd.concat([outstreamlines, df_streamlines], ignore_index=True)
        else:  # WITH tau decomposition
            # loop only over the coorect var
            nsteps = len(ds_ebm[loop_dim].values)
            nvals = ds_ebm[loop_dim].values
            # TODO ... deal with timing later
            kwargs['ts_year'] = "Agh.. need to fix timing later!"
            kwargs['slice_of_year'] ="fix please!"
            kwargs['yr_slices_remaining'] = nsteps
            kwargs['total_steps'] = nsteps
            # ...................................
            # start loop
            for thisstep in range(nsteps):
                kwargs['this_step'] = thisstep + 1
                ds_in = ds_ebm.sel({loop_dim: ds_ebm[loop_dim][thisstep]})
                # solve for tau
                ds_in, df_streamlines = Isotopes.taubar_withDecomposition(ds_in, iso_coords=df_coords, **kwargs)        
                # add info to streamlines
                df_streamlines[loop_dim] = nvals[thisstep]  # (the other dim is added w/in the taubar fxn)

                # create output data
                # ... SAVE RESULT
                ds_out = ds_in.assign_coords({loop_dim: nvals[thisstep]})
                ds_out = ds_out.expand_dims(dim=loop_dim)
                if thisstep == 0:
                    outds = ds_out.copy()
                    outstreamlines = df_streamlines
                else:
                    outds = xr.merge([outds, ds_out])
                    outstreamlines = pd.concat([outstreamlines, df_streamlines], ignore_index=True)
                

        # return result
        return outds, outstreamlines





    # [4] FULL LOOP
    @staticmethod
    def runloop(ds_clim_all, ds_force_all, df_coords=None, 
                runtype_list = ("full_model", "EBM_noIsotopes", "isotopes_noEBM"), **kwargs):

        # START THE RUN ----------------------------------
        l0 = "********************************************"
        l0a = ""
        l1 = "Now starting run " + kwargs['run_name']
        l1a = ""
        l2 = "--------------------------------------------"
        print(l0, l0a, l1, l1a, l2, sep=os.linesep)
        # ------------------------------------------------


        # set number of iterations
        slice_steps = kwargs.get('slices_to_solve')
        if ds_force_all is None:
            if 'time' in ds_clim_all.variables:
                timevalues = ds_clim_all.time.values
                tsteps = range(len(timevalues))
            else:
                timevalues = [1.]
                tsteps = range(1)
        else:
            timevalues = ds_force_all.time.values
            tsteps = range(len(timevalues))
        # if slice steps not defined, solve all
        if not slice_steps in ds_clim_all.yr_slice.values:
            slice_steps = ds_clim_all.yr_slice.values

        # check how to do any decomposition of tau_bar effects
        if kwargs.get('Tau-bar_decomp_local_v_regional_LocalEvapEffect') or kwargs.get('Tau-bar_decomp_local_v_regional_UpwindEffect'):
            spatial_decomposition = True
        else:
            spatial_decomposition = False
        clim_decomposition = kwargs.get('Tau-bar_decomp_E_L_s')
        tau_decompByTstep = kwargs.get('ELs_initstate_SameYrSlice')
        if tau_decompByTstep is None:
            tau_decompByTstep = 'no decomposition'

        
        # initialize time tracking
        kwargs['yr_slices_remaining_int'] = len(slice_steps) * len(timevalues) + 1  # add one b/c we subtract at start of tstep to keep all trackers together
        kwargs['yr_slices_remaining'] = str(kwargs.get('yr_slices_remaining_int'))
        

        # run the loop ------------------------

        outds, outstreamlines = Run.loop_ebm_iso(ds_clim_all, ds_force_all, slice_steps, tsteps, 
                                                 spatial_decomposition, clim_decomposition, tau_decompByTstep, 
                                                 df_coords=df_coords, **kwargs)
        
            
        # FINISH THE RUN ----------------------------------
        l0 = ""
        l1 = "********************************************"
        l2 = ""
        l3 = "run complete!                           :)"
        l4 = ""
        print(l0, l1, l2, l3, l4, l1, l0, sep=os.linesep)
        # ------------------------------------------------
        
        # return result
        return outds, outstreamlines







