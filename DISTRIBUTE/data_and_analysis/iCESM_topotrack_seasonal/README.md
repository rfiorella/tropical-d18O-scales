# README

This directory contains `*.nc` and `*.csv` files with iCESM results that include the isotope gradient data and relevant streamlines. 

1. `Climo_MonthlyMean_ALL.nc` includes the climatological monthly mean data for each relevant variable (unprocessed – it doesn't include upwind gradient data). 
2. `slice_*_ds_withDd18.nc` includes the climatology and spatial isotope gradient data for the relevant timeslices `ann_prWtd` (precipitation weighted annual mean), `djf` (december, january, february) and `jja` (june, july, august). 
3. `slice_*_streamlines.csv` same breakdown above, but with upwind streamline data. 
4. `slice_sznAnom_jja-djf_upwindCrossover.nc` includes the seasonal anomaly with the computed upwind crossover point. 


## Missing files (!!)
NOTE – this directory had 4 files larger than github's limit of 100MB. Please find those files separately on Zenodo: [Additional files for `spatial scales of tropical isotope change` project](https://zenodo.org/records/10836049?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjRiNTU5YzEwLWIxM2EtNDYwNC04YjBkLTIwY2ZlY2IxZmI3NSIsImRhdGEiOnt9LCJyYW5kb20iOiIzODM1MTFmZWJhMjQ3NzUxZWJhMGIyNjdiOTUxNTBhMiJ9.OEDPKuLeWStR8GBxXSJJl17aeD-dQMzoCUxKwgIskMY8icA6OcQXn-v3K36hP3fi34FK9AKINWMGMMKn569g5A). Contact me (tkukla13@gmail.org) with any questions or issues.