# README

This directory contains the data and processing scripts for Figure 3. It has three sub-directories with data, the processing script to produce them, and an accessory `*.nc` file containing data on the spatial isotope gradient anomaly. 

## `icesmdat` directory
Contains three sets of two files. In each set, `*_ae` means data are filtered by the "amount effect" (the precipitation anomaly is opposite in sign to the oxygen isotope anomaly):
1. `d18_anom*` is the seasonal isotope anomaly (y-axis in figure 3a) (per mille).
2. `d18Cross*` is the upwind distance where the d18O values intersect (x-axis in figure 3a) ("crossover") in the anomaly (km).
3. `Dd18_anom*` is the seasonal spatial isotope gradient anomaly (dot color in figure 3a) (per mille per 1000 km).

## `paleodat` directory
Contains three files with anomalies in the 218+GreenSah experiment and a directory (`2xCO2`) where those three files are repeated for the doubling CO2 experiment. All differences are precipitation-weighted mean annual in the experimental case minus the same in the CTRL case.
1. `paleo_dist.npy` is the upwind distance that causes 66% of the isotope anomaly (x-axis in figure 3b) (km).
2. `paleo_Dslope.npy` is the spatial isotope gradient anomaly (dot color in figure 3b) (per mille per 1000 km).
3. `paleo_Dtau.npy` is the change in attenuation (linearly related to d18O, y-axis in figure 3b). 

## `rtmdat` directory
Just one file (`rtm_df.csv`) which contains the x, y data for figure 3c. 