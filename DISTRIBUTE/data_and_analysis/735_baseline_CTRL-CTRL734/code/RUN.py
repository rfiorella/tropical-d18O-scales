# %% 
# --------------------------------------------------------
# MODEL TEST
# 
# --------------------------------------------------------

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

# MODEL 
from attenuationMod_fxns import *
from initialize import *

# %%
# read in input files -----------------------------
ds_clim_all, ds_force_all, df_coords, init = Run.inputfiles(**init)


# %%
# loop over time ----------------------------- 

outds, outstreamlines = Run.runloop(ds_clim_all, ds_force_all, df_coords, **init)


# ... save results ---------------------------------------------------------
results_dir = os.path.join(init['run_path'], init['run_name'], "results")
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
# name output files 
ds_loc = os.path.join(results_dir, init['run_name'] + "_CLIM.nc")
strmln_loc = os.path.join(results_dir, init['run_name'] + "_STREAMLINES.csv")
# save
outds.to_netcdf(ds_loc)
outstreamlines.to_csv(strmln_loc, index=False)
# --------------------------------------------------------------------------





