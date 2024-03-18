# README

This directory contains data and code to analyze the d18O response to orbital precession in the tropical SISALv2 database. 

## SISALv2 files 
1. `d13C.csv`
2. `d18o.csv`
3. `entity.csv`
4. `original_chronology.csv`
5. `sample.csv`
6. `site.csv`

## Analysis and results
1. `SISAL_smoothRange_perOrbCycle.R` is an R script that loops through each record, smooths it, and computes the d18O range for each orbital cycle. It outputs the results as the second file in this category.
2. `SISALv2_2kyrSmoothRangesINSOL_perOrbCycle.RDS` the results from running the first script. These data are used to produce the histogram in figure 1. 