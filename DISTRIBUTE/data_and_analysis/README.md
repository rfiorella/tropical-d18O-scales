# README

## Attenuation model runs
1. `735_baseline_CTRL-CTRL734` – the 2xCO2 simulation.
2. `733hiRes_218ka_367ppm_darkSah_baseline` – the 218+GreenSah simulation
3. `733a_hiRes_218ka_367ppm_darkSah_baseline_1000kmLocal` – the 218+GreenSah simulation repeated with a different distance (1000 km) for the local - non-local effect partitioning. 

## Other processing and data
1. `anom_v_distance` – contains data and code for figure 3; comparing the anomaly to its relevant upwind distance metric across three models (iCESM, attenuation model, rtm).
2. `iCESM_topotrack_seasonal` – contains the climatological monthly mean iCESM data, plus the processed timeslice files with upwind metrics like the isotope gradient and upwind topographic interaction. 
3. `rtm` – contains model scripts and code for running the reactive transport model.
4. `sisal_v2_d18Oranges` – contains data and scripts for analyzing the SISALv2 tropical data on precession timescales. 