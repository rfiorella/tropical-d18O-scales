# -------------------------------------------------
# 
# accessory functions that are fed into the 
#   main model code 
# 
# -------------------------------------------------

import numpy as np
import xarray as xr
import pandas as pd
import math
from scipy.interpolate import RegularGridInterpolator  # for tau-bar computation
import random


# [0] BASIC FUNCTIONS
def var_check(vars, dataset):
    var_check_result = all(var in dataset.variables for var in vars)
    if not var_check_result:
        joined_strings = ", ".join(vars)
        raise ValueError(f"Missing or mis-named variables. Required vars: {joined_strings}")

# FUNCTION TO convert from grid points to lat lon (required for stitch contours fxn)
# ... idx_input is the 1d array of values you want to translate to degrees
# ... deg_grid is the ordered 1d array where index 0 equals min(deg_grid) and index of len(deg_grid) equals max(deg_grid) 
def linear_idx_to_degrees(idx_input, deg_grid):
    '''
    Converts a 1-d index [0,n] to values in degrees using y = mx + b

    Parameters
    ----------

    idx_input : 1d array
            indices denoting the ordered lat or lon 
            (0 must be minimum lat or lon, and len(deg_grid) must be max)

    deg_grid : 1d array
            ordered where 0 equals min(deg_grid) and index of len(deg_grid) 
            equals max(deg_grid)
    '''


    degmin = deg_grid.min() ; degmax = deg_grid.max()
    idxmin = 0 ; idxmax = len(deg_grid)
    # slope and itnercept
    m_slope = (degmax - degmin) / (idxmax - idxmin)  # slope of linear relationship
    b_int = degmax - (m_slope * idxmax)
    # y=mx + b
    out_y = m_slope * idx_input + b_int
    return out_y



# [1] FUNCTIONS WITHIN TAUBAR
# tau and climatology streamlines for a given point 
def tau_streamline_point(lat0, lon0, Nmax, taumax, Dx, dx, Efit, Pfit, 
                         vqfit, uqfit, Fmagfit, lfracfit):
    # initialize arrays
    tau = np.zeros(Nmax)
    E0 = np.zeros(Nmax)
    P0 = np.zeros(Nmax)
    PminE = np.zeros(Nmax)  # diagnostic - not used in calculation
    Fmag0 = np.zeros(Nmax)
    dist = np.zeros(Nmax)   # for tracking 1d evolution (km)
    vq0 = np.zeros(Nmax)
    uq0 = np.zeros(Nmax)
    lfrac0 = np.zeros(Nmax)
    
    # array for spatial evolution of tau-bar for diagnostics
    # tau_bar_streamline = np.zeros(Nmax)  # for spatial evolution of tau-bar, calculated on
    # check_last_km = 20e3   # [km] tau_bar calculated only over this downwind region (to avoid issues at "source")
    # get_taubar_streamline = False  # if true, run a loop to get tau-bar for end of streamline

    # INTERPOLATE to point
    E0[0] = Efit((lat0,lon0)) # interpolate evaporation
    P0[0] = Pfit((lat0,lon0)) # interpolate precipitation
    PminE[0] = P0[0] - E0[0]  # diagnostic - not used in calculation
    vq0[0] = vqfit((lat0,lon0))
    uq0[0] = uqfit((lat0,lon0))
    Fmag0[0] = Fmagfit((lat0,lon0)) # interpolate vapor transport
    lfrac0[0] = lfracfit((lat0, lon0))
    dist[0] = 0   # start at zero km
    dist_inland = 0  # this gets updated throughout the loop
    coast_step = 0   # when we find the coast, we note the step-number for tau-adjuster 

    # save lat/lon path for troubleshooting
    lat_save = np.zeros(Nmax) ; lon_save = np.zeros(Nmax)
    lat_save[0] = lat0 ; lon_save[0] = lon0

    # MARCH UPSTREAM
    # march upstream and find all variables along vapor transport path
    # end when tau >= taumax or when distance >= dmax.
    step = 0
    this_landmass = True  # track if we're still over the same landmass as the sink
    while tau[step] < taumax and step < Nmax-1: # should change to Nmax-1 but then it takes forever??
        
        # find incremental change in lat/lon based on uq and vq (in degrees)
        #dtheta = np.logaddexp(-vq0,-Fmag0[i]*Dx)
        dtheta = -vq0[step]/Fmag0[step]*Dx 
        #dphi = np.logaddexp(-uq0,-Fmag0[i]*Dx)/np.cos(np.radians((lat0+dtheta/2)))
        # avoid edge effects at singularity of 90 deg
        if lat0+dtheta > 89.5:
            dtheta = 89.5 - lat0
        if lat0+dtheta < -89.5:
            dtheta = -89.5 - lat0
        dphi = -uq0[step]/Fmag0[step]*Dx/np.cos(np.radians((lat0+dtheta/2)))

        # get distance in this step
        xdist = dx ; ydist = dx/np.cos(np.radians((lat0+dtheta/2)))
        dist[step+1] = dist[step] + np.sqrt(xdist**2 + ydist**2)  # distance between steps

        lat1 = lat0+dtheta
        lon1 = lon0+dphi

        # ensure that latitude stays within -90<=lat<=90
        if lat1 < -90:
            lat1 = 180+lat1
            lon1 = lon1-180
        elif lat1 > 90:
            lat1 = 180-lat1
            lon1 = lon1-180
        
        # ensure that longitude stays within 0<=lon<=360
        if lon1 < 0:
            lon1 = lon1+360
        elif lon1 > 360:
            lon1 = lon1-360

        # save for troubleshooting
        lat_save[step+1] = lat1 ; lon_save[step+1] = lon1
        
        vq0[step+1] = vqfit((lat1,lon1))
        uq0[step+1] = uqfit((lat1,lon1))
        E0[step+1] = Efit((lat1,lon1))
        Fmag0[step+1] = Fmagfit((lat1,lon1))
        P0[step+1] = Pfit((lat1,lon1))
        lfrac0[step+1] = lfracfit((lat1,lon1))
        PminE[step+1] = P0[step+1] - E0[step+1]  # diagnostic (not used in calcs)

        # add to distance if lfrac
        if lfrac0[step+1] > 0.6 and this_landmass:
            dist_inland = dist[step+1]
        if lfrac0[step+1] <= 0.6:  # if we've stepped upstream over ocean, don't count dist inland any more
            this_landmass = False
            if coast_step == 0:    # if we haven't noted the coast idx, note it now (then not again)
                coast_step = step

        # trapezoidal integration of mu; multiply dx by 1000 to convert to meters
        tau[step+1] = tau[step]+(dx*1000)*(P0[step]/Fmag0[step]+P0[step+1]/Fmag0[step+1])/2

        # vq0 = vq1
        # uq0 = uq1
        lat0 = lat1
        lon0 = lon1

            
        step = step+1 # step forward
    
    # handle special case where air-mass never goes over ocean
    if dist_inland > 0 and coast_step == 0:
        coast_step = step   # set the "coast" to the last step, and we just integrate the whole streamline

    # where it's not nan
    nanindex = tau > 0; nanindex[0] = 1 

    # interpolate any E0 nans (not sure why they sometimes appear)
    e0nans, e02 = np.isnan(E0), lambda z: z.nonzero()[0]
    if len(np.argwhere(e0nans)) > 0:  # interpolate to remove nans
        E0[e0nans] = np.interp(e02(e0nans), e02(~e0nans), E0[~e0nans])

    # additional terms
    mu = P0[nanindex] / Fmag0[nanindex]
    wp = (E0[nanindex] * math.e**(-tau[nanindex])) / (np.trapz(E0[nanindex] * math.e**(-tau[nanindex])))

    # return results
    return coast_step, nanindex, dist, dist_inland, tau, mu, wp, E0, P0, vq0, uq0, Fmag0, lfrac0, PminE, lat_save, lon_save
    

# initialize arrays for interpolation
def tau_arrayInitialize(ds, lat, lon, **kwargs):
    # collect names
    decomposition_dim = kwargs.get("decomp")
    U_name = kwargs.get("ZonalQFlux_Field")
    V_name = kwargs.get("MeridQFlux_Field")
    e_name = kwargs.get("Evaporation_Field")
    p_name = kwargs.get("Precipitation_Field")
    fmag_name = "Fmag"
    if e_name is None:
        e_name = "ET"
    if p_name is None:
        p_name = "PRECT"

    if decomposition_dim is None:  # only build lat / lon arrays
        # create empty variables
        P = np.zeros((len(lat),len(lon)))
        E = np.zeros_like(P) ; Fmag = np.zeros_like(P)
        UQ = np.zeros_like(P) ; VQ = np.zeros_like(P)
        lfrac = np.zeros_like(P)

        # pull out variables
        E[:,:], P[:,:] = ds[e_name], ds[p_name]
        UQ[:,:], VQ[:,:] = ds[U_name], ds[V_name]
        Fmag[:,:] = ds[fmag_name]
        lfrac[:,:] = ds['LANDFRAC']
    else:
        # create empty variables
        P = np.zeros((len(lat),len(lon), len(decomposition_dim)))
        E = np.zeros_like(P) ; Fmag = np.zeros_like(P)
        UQ = np.zeros_like(P) ; VQ = np.zeros_like(P)
        lfrac = np.zeros_like(P)

        # pull out variables
        # (depending on the decomposition_dim, they might need to be transposed)
        E[:,:,:] = ds[e_name].transpose('lat', 'lon', decomposition_dim.name)
        P[:,:,:] = ds[p_name].transpose('lat', 'lon', decomposition_dim.name)
        UQ[:,:,:] = ds[U_name].transpose('lat', 'lon', decomposition_dim.name)
        VQ[:,:,:] = ds[V_name].transpose('lat', 'lon', decomposition_dim.name)
        Fmag[:,:,:] = ds[fmag_name].transpose('lat', 'lon', decomposition_dim.name)
        lfrac[:,:,:] = ds['LANDFRAC'].transpose('lat', 'lon', decomposition_dim.name)

    # return arrs
    return P, E, UQ, VQ, Fmag, lfrac



# build tau grid inteprolators for the streamline computation
def tau_interpolatorInitialize(lon, lat, E, P, UQ, VQ, lfrac):
    # first check if lat is in xr or np arr
    # (get into array if it's not already)
    if not isinstance(lat, np.ndarray):
        lat = lat.values
    if not isinstance(lon, np.ndarray):
        lon = lon.values
    # create buffer so that data spans all longitude points for interpolation
    #lat_test = np.insert(lat,0,lat[-1]-180); lat_test = np.append(lat_test,lat[0]+180)
    lon_test = np.insert(lon,0,lon[-1]-360); lon_test = np.append(lon_test,lon_test[1]+360)
    E_test = np.concatenate((E[:,-1,None],E, E[:,0,None]),1)
    P_test = np.concatenate((P[:,-1,None],P,P[:,0,None]),1)
    UQ_test = np.concatenate((UQ[:,-1,None],UQ,UQ[:,0,None]),1)
    VQ_test = np.concatenate((VQ[:,-1,None],VQ,VQ[:,0,None]),1)
    Fmag_test = np.sqrt(UQ_test**2+VQ_test**2) # find magnitude of transport vectors
    lfrac_test = np.concatenate((lfrac[:,-1,None], lfrac, lfrac[:,0,None]),1)

    # create interpolation functions for all variables
    Efit = RegularGridInterpolator((lat,lon_test),E_test,method='linear')
    Pfit = RegularGridInterpolator((lat,lon_test),P_test,method='linear')
    uqfit = RegularGridInterpolator((lat,lon_test),UQ_test,method='linear')
    vqfit = RegularGridInterpolator((lat,lon_test),VQ_test,method='linear')
    Fmagfit = RegularGridInterpolator((lat,lon_test),Fmag_test,method='linear')
    lfracfit = RegularGridInterpolator((lat, lon_test), lfrac_test, method ='linear')

    # return results
    return Efit, Pfit, uqfit, vqfit, Fmagfit, lfracfit



# add data to a streamlines dataframe
def build_streamlines_df(lat, lon, x, y, idx_counter, no_streamlines_yet, streamline_collect, 
                         save_streamline, dist, lat_save, lon_save, E0, Fmag0, P0, PminE, 
                         tau, mu, wp, lfrac0, nanindex, streamline_type,
                         df_streamlines=None, **kwargs):
    # first check if lat is in xr or np arr
    # (get into array if it's not already)
    if not isinstance(lat, np.ndarray):
        lat = lat.values
    if not isinstance(lon, np.ndarray):
        lon = lon.values
    # figure out what we need to ignore 
    ignore_these = kwargs.get('streamline_ignore_vars')
    if ignore_these is not None:
        if "nanindex" in ignore_these:
            nanindex = np.repeat(True, len(tau))
        if "lat_save" in ignore_these:
            lat_save = np.repeat(-9999, len(tau))
        if "lon_save" in ignore_these:
            lon_save = np.repeat(-9999, len(tau))
        if "lfrac0" in ignore_these:
            lfrac0 = np.repeat(-9999, len(tau))
    # output that we need no matter what (especially if we compute tau, but don't save streamline and we want to decompose tau)
    tau_out = tau[nanindex]
    E0_out = E0[nanindex]
    dist_out = dist[nanindex]
    mu_out = mu

    # now we build streamlines
    if streamline_collect in ['all', 'some']:
        if save_streamline:
            dist_x = dist[nanindex]
            lat_save_x = lat_save[nanindex]
            lon_save_x = lon_save[nanindex]
            E0_x = E0[nanindex]
            Fmag_x = Fmag0[nanindex]
            P0_x = P0[nanindex]
            PminE_x = PminE[nanindex]
            tau_x = tau[nanindex]
            mu_x = mu    # nanindex applied earlier
            wp_x = wp    # nanindex applied earlier
            lfrac_x = lfrac0[nanindex]
            # define nx
            nx = kwargs.get('StreamlineSave_Coarsener')
            if nx is None:
                nx = 1
        # set to na value otherwise
        else: # (only one of the -9999s will be saved)
            dist_x = np.array([-9999, -9999])
            lat_save_x = np.array([-9999, -9999])
            lon_save_x = np.array([-9999, -9999])
            E0_x = np.array([-9999, -9999])
            Fmag_x = np.array([-9999, -9999])
            P0_x = np.array([-9999, -9999])
            PminE_x = np.array([-9999, -9999])
            tau_x = np.array([-9999, -9999])
            mu_x = np.array([-9999, -9999])
            wp_x = np.array([-9999, -9999])
            lfrac_x = np.array([-9999, -9999])
            # change nx
            nx = 2

        if no_streamlines_yet:
            df_streamlines = {
                "dist_km": dist_x[::nx],
                "lat_streamline": lat_save_x[::nx],
                "lon_streamline": lon_save_x[::nx],
                "ET": E0_x[::nx],
                "Fmag": Fmag_x[::nx],
                "PRECT": P0_x[::nx],
                "PminE": PminE_x[::nx],
                "tau": tau_x[::nx],
                "mu": mu_x[::nx],
                "wp": wp_x[::nx],
                "lfrac": lfrac_x[::nx]
            }
            df_streamlines['lat_sink'] = lat[y]
            df_streamlines['lon_sink'] = lon[x]
            df_streamlines['streamline_type'] = streamline_type
            df_streamlines['idx'] = idx_counter
            df_streamlines = pd.DataFrame(df_streamlines)
        else:
            tmpdf_streamlines = {
                "dist_km": dist_x[::nx],
                "lat_streamline": lat_save_x[::nx],
                "lon_streamline": lon_save_x[::nx],
                "ET": E0_x[::nx],
                "Fmag": Fmag_x[::nx],
                "PRECT": P0_x[::nx],
                "PminE": PminE_x[::nx],
                "tau": tau_x[::nx],
                "mu": mu_x[::nx],
                "wp": wp_x[::nx],
                "lfrac": lfrac_x[::nx]
            }
            tmpdf_streamlines['lat_sink'] = lat[y]
            tmpdf_streamlines['lon_sink'] = lon[x]
            df_streamlines['streamline_type'] = streamline_type
            tmpdf_streamlines['idx'] = idx_counter
            tmpdf_streamlines = pd.DataFrame(tmpdf_streamlines)
            # combine
            df_streamlines = pd.concat([df_streamlines, tmpdf_streamlines], ignore_index=True)
    else:
        df_streamlines = {
                "dist_km": 0,
                "lat_streamline": 0,
                "lon_streamline": 0,
                "ET": 0,
                "Fmag": 0,
                "PRECT": 0,
                "PminE": 0,
                "tau": 0,
                "mu": 0,
                "wp": 0,
                "lfrac": 0
        }
        df_streamlines['lat_sink'] = 0
        df_streamlines['lon_sink'] = 0
        df_streamlines['streamline_type'] = streamline_type
        df_streamlines['idx'] = 0
        df_streamlines = pd.DataFrame(df_streamlines)

    return df_streamlines, tau_out, E0_out, dist_out, mu_out


# cast in tau coordinates
def dist_to_tau_coords(tau, nanindex, mu, wp, E0):
    # (this is necessary to decompose E and L)
    tau_max, tau_min = np.max(tau[nanindex]), np.min(tau[nanindex])
    n_tausteps = len(tau[nanindex])
    tausteps = np.linspace(tau_min, tau_max, n_tausteps)
    # linear interpolation
    mu_taugrid = np.interp(tausteps, tau[nanindex], mu)
    wp_taugrid = np.interp(tausteps, tau[nanindex], wp)
    E0_taugrid = np.interp(tausteps, tau[nanindex], E0[nanindex])
    # return results
    return mu_taugrid, wp_taugrid, E0_taugrid, tausteps


# decompose local v regional tau effects
def tau_spatial_decompose(local_threshold, x, y, thisidx,
                          dist, tau, E0, nanindex,
                          disti, taui, E0i, nanindexi,
                          Dtau_bar_regionalEff=None, Dtau_bar_localEff=None,
                          localEffMode = False):
    # local and remote values for tau computation
    local_idx_old = np.where(disti[nanindexi] <= local_threshold)[0]
    local_idx_new = np.where(dist[nanindex] <= local_threshold)[0]
    # naninex maxis
    nandxmax = np.nanmax(np.where(nanindex)[0])
    nandxmaxi = np.nanmax(np.where(nanindexi)[0])
    # old local v regional
    tau_old_local, tau_old_regional = taui[local_idx_old], taui[np.max(local_idx_old)+1:nandxmaxi]
    E0_old_local, E0_old_regional = E0i[local_idx_old], E0i[np.max(local_idx_old)+1:nandxmaxi]
    # new local v regional
    tau_new_local, tau_new_regional = tau[local_idx_new], tau[np.max(local_idx_new)+1:nandxmax]
    E0_new_local, E0_new_regional = E0[local_idx_new], E0[np.max(local_idx_new)+1:nandxmax]
    # get local vs regional constant inputs for tau computation
    tau_localConst = np.concatenate((tau_old_local, tau_new_regional))
    tau_regionalConst = np.concatenate((tau_new_local, tau_old_regional))
    E0_localConst = np.concatenate((E0_old_local, E0_new_regional))
    E0_regionalConst = np.concatenate((E0_new_local, E0_old_regional))
    # local constant tau_bar
    taubar_localConst = np.trapz(tau_localConst*E0_localConst*math.e**(-tau_localConst))/np.trapz(E0_localConst*math.e**(-tau_localConst))
    # regional constant tau_bar
    taubar_regionalConst = np.trapz(tau_regionalConst*E0_regionalConst*math.e**(-tau_regionalConst))/np.trapz(E0_regionalConst*math.e**(-tau_regionalConst))
    # get initial tau
    taubar_old = np.trapz(taui[nanindexi]*E0i[nanindexi]*math.e**(-taui[nanindexi]))/np.trapz(E0i[nanindexi]*math.e**(-taui[nanindexi]))
    
    # local effect mode is used to return just one value (for use with dist_to_Dtaubar_frac function)
    if localEffMode:
        localEff = taubar_regionalConst - taubar_old
        return localEff
    else:
        # and the change due to local effects
        Dtau_bar_regionalEff[y,x,thisidx] = taubar_localConst - taubar_old 
        Dtau_bar_localEff[y,x,thisidx] = taubar_regionalConst - taubar_old

        return Dtau_bar_regionalEff, Dtau_bar_localEff


# get distance from sink to the point where some fraction of Dtau_bar is accounted for
# ... this function computes tau old and tau new and splices them together, which
# ... means that changes in local P do not affect the isotopic composition of remotely
# ... sourced (evaporated) moisture. Basically, this function analyzes the effect of 
# ... changes in local P on local moisture source
#     (in the next fxn, we get old and new mu, and get new tau-bar from a new tau for
#      each case, such that changes in local P act on remotely sourced E)
def dist_to_Dtaubar_frac_LocalEvap(x, y, thisidx, Dtau_bar_50_dist, Dtau_bar, 
                         dist, tau, E0, nanindex,
                         disti, taui, E0i, nanindexi,
                         Dtau_frac, Dtau_bar_dist_percError = 1, troubleshoot=False):
    this_Dtau = Dtau_bar[y,x,thisidx]
    this_Dtau_frac = this_Dtau * Dtau_frac
    # get a range of allowable Dtau values (when we find a distance within this range, we stop)
    # ... if Dtau_bar_dist_percError is 5, then we will take any fractional change within 5% of the target 
    Dtau_frac_range = [this_Dtau_frac - (Dtau_bar_dist_percError/100) * this_Dtau_frac, this_Dtau_frac + (Dtau_bar_dist_percError/100) * this_Dtau_frac]
    min_Dtau_frac = np.min(Dtau_frac_range)
    max_Dtau_frac = np.max(Dtau_frac_range)
    # get distance range 
    mindist, maxdist = 0, int(np.nanmax(dist))
    # initialize while loop
    if troubleshoot:  # then save arrays, rather than update the one value
        localEff_arr=np.zeros(50000) ; localEff_arr[0] = -1e3
        testdistance_arr = np.zeros(50000) ; testdistance_arr[0] = -1e3
    localEff = -1e3   # outrageous starting point (just any number that will never reasonably be w/in min and max Dtau_frac)
    testdistance = 0 # negative so we know it won't work (it gets updated right away)
    searchiter = 0
    binary_search_mode = True  # if we can't find a solution, we do a brute force search before quitting
    brute_force_search_step = 5  # step size for brute force search (most grid cells don't require this)
    searchiterMax = 500
    # while loop -- while localEff is not within the accepted effect range, keep looking
    # ... performs binary search
    while localEff < min_Dtau_frac or localEff > max_Dtau_frac:
        # test distance splits min and max
        # testdistance += 10 # (troubleshoot)
        if searchiter == 0:
            testdistance = 1
        elif binary_search_mode: 
            testdistance = (mindist + maxdist) // 2
        else:
            testdistance += brute_force_search_step

        # run spatial decompose fxn
        localEff = tau_spatial_decompose(testdistance, x, y, thisidx,
                                        dist, tau, E0, nanindex,
                                        disti, taui, E0i, nanindexi,
                                        localEffMode = True)
        if troubleshoot:
            testdistance_arr[searchiter] = testdistance
            localEff_arr[searchiter] = localEff
        # update min and max
        # if not far enough from zero, then look further away
        if np.abs(localEff) < np.abs(min_Dtau_frac): # then increase the distance
            mindist = testdistance  # update minimum allowable dist
        # if too far, then look closer 
        if np.abs(localEff) > np.abs(max_Dtau_frac): # then decrease the distance
            maxdist = testdistance  # update maximum allowable dist
        
        # update search iter
        searchiter += 1
        # break the loop if it goes too long
        if searchiter == searchiterMax and binary_search_mode: # arbitrarily large number
            testdistance=10 # reset test distance
            # reset 
            if troubleshoot:
                localEff_arr=np.zeros_like(localEff_arr) ; localEff_arr[0] = -1e3
                testdistance_arr = np.zeros_like(testdistance_arr) ; testdistance_arr[0] = -1e3
            searchiter = 0
            searchiterMax = int(np.nanmax(dist)/brute_force_search_step)
            binary_search_mode = False
        elif searchiter >= searchiterMax:
            print("Couldn't find a distance solution in 'fxn.dist_to_Dtaubar_frac' ")
            testdistance = np.nan
            break

    # assign the final answer
    Dtau_bar_50_dist[y,x,thisidx] = testdistance
    # return result
    return Dtau_bar_50_dist



# create composite tau streamline from local and remote effects (if tau_local_v_regional and UpwindEff are True)
def local_plus_regional_tauStreamline(local_threshold, dx,
                                        P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
                                        P0, Fmag0, E0, dist, nanindex, lfrac0):
    # get composite local and remote streamlines
    local_idx_old = np.where(disti[nanindexi] <= local_threshold)[0]
    local_idx_new = np.where(dist[nanindex] <= local_threshold)[0]
    # naninex maxis
    nandxmax = np.nanmax(np.where(nanindex)[0])
    nandxmaxi = np.nanmax(np.where(nanindexi)[0])
    # old local v regional
    P0_old_local, P0_old_regional = P0i[local_idx_old], P0i[np.max(local_idx_old)+1:nandxmaxi]
    Fmag0_old_local, Fmag0_old_regional = Fmag0i[local_idx_old], Fmag0i[np.max(local_idx_old)+1:nandxmaxi]
    E0_old_local, E0_old_regional = E0i[local_idx_old], E0i[np.max(local_idx_old)+1:nandxmaxi]
    # new local v regional
    P0_new_local, P0_new_regional = P0[local_idx_new], P0[np.max(local_idx_new)+1:nandxmax]
    Fmag0_new_local, Fmag0_new_regional = Fmag0[local_idx_new], Fmag0[np.max(local_idx_new)+1:nandxmax]
    E0_new_local, E0_new_regional = E0[local_idx_new], E0[np.max(local_idx_new)+1:nandxmax]
    # distance - new
    dist_old_local, dist_old_regional = disti[local_idx_old], disti[np.max(local_idx_old)+1:nandxmaxi]
    dist_new_local, dist_new_regional = dist[local_idx_new], dist[np.max(local_idx_new)+1:nandxmax]

    # get local vs regional constant inputs for tau computation
    P0_localConst = np.concatenate((P0_old_local, P0_new_regional))
    P0_regionalConst = np.concatenate((P0_new_local, P0_old_regional))
    Fmag0_localConst = np.concatenate((Fmag0_old_local, Fmag0_new_regional))
    Fmag0_regionalConst = np.concatenate((Fmag0_new_local, Fmag0_old_regional))
    E0_localConst = np.concatenate((E0_old_local, E0_new_regional))
    E0_regionalConst = np.concatenate((E0_new_local, E0_old_regional))
    dist_localConst = np.concatenate((dist_old_local, dist_new_regional))
    dist_regionalConst = np.concatenate((dist_new_local, dist_old_regional))

    # integrate each step
    tau_localConst = np.zeros(len(P0_localConst))
    tau_regionalConst = np.zeros(len(P0_regionalConst))
    max_steps = np.max((len(P0_regionalConst), len(P0_localConst)))
    for step in range(max_steps):
        if step < (len(P0_localConst)-1): # then we can solve for tau
            # trapezoidal integration of mu; multiply dx by 1000 to convert to meters
            tau_localConst[step+1] = tau_localConst[step]+(dx*1000)*(P0_localConst[step]/Fmag0_localConst[step]+P0_localConst[step+1]/Fmag0_localConst[step+1])/2
        if step < (len(P0_regionalConst)-1): # then we can solve for tau
            # trapezoidal integration of mu; multiply dx by 1000 to convert to meters
            tau_regionalConst[step+1] = tau_regionalConst[step]+(dx*1000)*(P0_regionalConst[step]/Fmag0_regionalConst[step]+P0_regionalConst[step+1]/Fmag0_regionalConst[step+1])/2

    # return results
    return tau_localConst, E0_localConst, tau_regionalConst, E0_regionalConst, Fmag0_localConst, Fmag0_regionalConst, P0_localConst, P0_regionalConst, dist_localConst, dist_regionalConst

    
# distance to some fraction of the isotope anomaly being explained when 
# accounting for the streamline / upwind effect of more rainout 
def dist_to_Dtaubar_frac_streamline(x, y, thisidx, dx, Dtau_bar, taubar_init,
                                        P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
                                        P0, Fmag0, E0, dist, nanindex, lfrac0,
                                        Dtau_bar_50_dist_strmline,
                                Dtau_frac, Dtau_bar_dist_percError = 1, troubleshoot=False):
    this_Dtau = Dtau_bar[y,x,thisidx]
    this_Dtau_frac = this_Dtau * Dtau_frac
    # get a range of allowable Dtau values (when we find a distance within this range, we stop)
    # ... if Dtau_bar_dist_percError is 5, then we will take any fractional change within 5% of the target 
    Dtau_frac_range = [this_Dtau_frac - (Dtau_bar_dist_percError/100) * this_Dtau_frac, this_Dtau_frac + (Dtau_bar_dist_percError/100) * this_Dtau_frac]
    min_Dtau_frac = np.min(Dtau_frac_range)
    max_Dtau_frac = np.max(Dtau_frac_range)
    # get distance range 
    mindist, maxdist = 0, int(np.nanmax(dist))
    # initialize while loop
    if troubleshoot:  # then save arrays, rather than update the one value
        localEff_arr=np.zeros(50000) ; localEff_arr[0] = -1e3
        testdistance_arr = np.zeros(50000) ; testdistance_arr[0] = -1e3
    localEff = -1e3   # outrageous starting point (just any number that will never reasonably be w/in min and max Dtau_frac)
    testdistance = 0 # negative so we know it won't work (it gets updated right away)
    searchiter = 0
    binary_search_mode = True  # if we can't find a solution, we do a brute force search before quitting
    brute_force_search_step = 5  # step size for brute force search (most grid cells don't require this)
    searchiterMax = 500
    
    # while loop -- while localEff is not within the accepted effect range, keep looking
    # ... performs binary search
    while localEff < min_Dtau_frac or localEff > max_Dtau_frac:
        # test distance splits min and max
        # testdistance += 10 # (troubleshoot)
        if searchiter == 0:
            testdistance = 1
        elif binary_search_mode: 
            testdistance = (mindist + maxdist) // 2
        else:
            testdistance += brute_force_search_step

        # get the local / regional effects
        tau_localConst, E0_localConst, tau_regionalConst, E0_regionalConst, _, _, _, _, _, _  = local_plus_regional_tauStreamline(testdistance, dx,
                                                                                                                P0i, Fmag0i, E0i, disti, nanindexi, lfrac0i,
                                                                                                                P0, Fmag0, E0, dist, nanindex, lfrac0)
        # get tau_bar due to regional change
        taubar_localConst = np.trapz(tau_localConst*E0_localConst*math.e**(-tau_localConst))/np.trapz(E0_localConst*math.e**(-tau_localConst))
        # and due to local change
        taubar_regionalConst = np.trapz(tau_regionalConst*E0_regionalConst*math.e**(-tau_regionalConst))/np.trapz(E0_regionalConst*math.e**(-tau_regionalConst))
        # get local effect
        localEff = taubar_regionalConst - taubar_init
        
        if troubleshoot:
            testdistance_arr[searchiter] = testdistance
            localEff_arr[searchiter] = localEff
        # update min and max
        # if not far enough from zero, then look further away
        if np.abs(localEff) < np.abs(min_Dtau_frac): # then increase the distance
            mindist = testdistance  # update minimum allowable dist
        # if too far, then look closer 
        if np.abs(localEff) > np.abs(max_Dtau_frac): # then decrease the distance
            maxdist = testdistance  # update maximum allowable dist
        
        # update search iter
        searchiter += 1
        # break the loop if it goes too long
        if searchiter == searchiterMax and binary_search_mode: # arbitrarily large number
            testdistance=10 # reset test distance
            # reset 
            if troubleshoot:
                localEff_arr=np.zeros_like(localEff_arr) ; localEff_arr[0] = -1e3
                testdistance_arr = np.zeros_like(testdistance_arr) ; testdistance_arr[0] = -1e3
            searchiter = 0
            searchiterMax = int(np.nanmax(dist)/brute_force_search_step)
            binary_search_mode = False
        elif searchiter >= searchiterMax:
            print("Couldn't find a distance solution in 'fxn.dist_to_Dtaubar_frac' ")
            testdistance = np.nan
            break
    
    # return the final answer
    Dtau_bar_50_dist_strmline[y,x,thisidx] = testdistance
    # return result
    return Dtau_bar_50_dist_strmline




# decompose E, L, path tau effects
def tau_decompose_ELPath(x, y, thisidx, taubar_init,
                         E0_E, E0_L, E0_s,
                         mu_E, mu_L, mu_s,
                         wp_E, wp_L, wp_s,
                         tau_E, tau_L, tau_s, 
                         nanindex_E, nanindex_L, nanindex_s,
                         Dtau_bar_newE, Dtau_bar_newL, Dtau_bar_newPath):
    # get the decomposed changes in taubar
    # Dtau -- new E
    mu_taugrid_E, _, E0_taugrid_E, tausteps_E = dist_to_tau_coords(tau_E, nanindex_E, mu_E, wp_E, E0_E)
    tau_bar_newE =  np.trapz(E0_taugrid_E * (1/mu_taugrid_E) * tausteps_E * math.e**(-tausteps_E)) / np.trapz(E0_taugrid_E * (1/mu_taugrid_E) * math.e**(-tausteps_E))
    Dtau_bar_newE[y,x,thisidx] = tau_bar_newE - taubar_init
    # Dtau -- new L
    mu_taugrid_L, _, E0_taugrid_L, tausteps_L = dist_to_tau_coords(tau_L, nanindex_L, mu_L, wp_L, E0_L)
    tau_bar_newL =  np.trapz(E0_taugrid_L * (1/mu_taugrid_L) * tausteps_L * math.e**(-tausteps_L)) / np.trapz(E0_taugrid_L * (1/mu_taugrid_L) * math.e**(-tausteps_L))
    Dtau_bar_newL[y,x,thisidx] = tau_bar_newL - taubar_init
    # Dtau -- new Path
    mu_taugrid_s, _, E0_taugrid_s, tausteps_s = dist_to_tau_coords(tau_s, nanindex_s, mu_s, wp_s, E0_s)
    tau_bar_news =  np.trapz(E0_taugrid_s * (1/mu_taugrid_s) * tausteps_s * math.e**(-tausteps_s)) / np.trapz(E0_taugrid_s * (1/mu_taugrid_s) * math.e**(-tausteps_s))
    Dtau_bar_newPath[y,x,thisidx] = tau_bar_news - taubar_init

    # return result
    return Dtau_bar_newE, Dtau_bar_newL, Dtau_bar_newPath



# compute the terrestrial fraction of evaporation source
def terrestrial_E_frac(land_check, wp, coast_step, tau, E0, nanindex, n_samples):
    # landfrac_streamline = dist_inland / np.max(dist[nanindex])
    E_terrFrac = np.sum(wp[:coast_step])
    if land_check > 0.6: # if this is terrestrial point, get the e-wtd tauba
        land_wp = wp[:coast_step]
        land_Esource_frac = np.trapz(land_wp) / np.trapz(wp) # get land_wp over total
        nsteps = len(land_wp)
        
        if n_samples > len(land_wp):  # don't sample things that don't exist
            n_samples = len(land_wp)
        # compute sampling step
        samples_idx = np.linspace(0, nsteps - 1, n_samples, dtype=int)
        # loop through all samples 
        this_wp = np.empty(n_samples)
        this_taubar = np.empty_like(this_wp)
        idx = 0
        for thissample in samples_idx:
            this_wp[idx] = land_wp[thissample]
            taunan = tau[nanindex] ; E0nan = E0[nanindex]
            this_taubar[idx] = np.trapz(taunan[thissample:]*E0nan[thissample:]*math.e**(-taunan[thissample:]))/np.trapz(E0nan[thissample:]*math.e**(-taunan[thissample:]))
            idx += 1
        tau_wtd = this_taubar * this_wp
        total_wp = np.sum(this_wp)
        tau_wtd_mean = np.sum(tau_wtd) / total_wp
    else:
        tau_wtd_mean = 0
        land_Esource_frac = 0

    return tau_wtd_mean, land_Esource_frac, E_terrFrac