# README

This directory contains [code](https://github.com/tykukla/Vapor_Transport_Model_KuklaEtAl2019) to run the reactive transport model [Kukla et al., 2019; JGR Atmos](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JD029571) plus the outputs used in the present study. 

## RTM scripts 
See [Kukla et al., 2019](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JD029571) for more information.
1. `build_topo_fxn.R`
2. `constants.R`
3. `hydro_fxns.R`
4. `isofrac_fxns.R`
5. `ModelRUN.R`
6. `ModelSolve_fxn.R`
7. `ODE_fxns.R`

## Analyses
These scripts run the two model analyses presented in the text by modifying the `ModelRUN.R` script:
1. `ModelRUN_UpwindMchange+Mchange_v2.R` runs a series of simulations varying the upwind rainout distance and the rainout intensity at the same time, producing results that are plotted in figure 5.
2. `ModelRUN_UpwindMchange.R` runs a series of simulations where only the upwind rainout distance is varied, producing results plotted in figure 3c. 

Model results for each are found in both `*.RDS` and `*.csv` form.