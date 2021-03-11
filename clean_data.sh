#!/usr/bin/bash

# Run each data cleaning script in its own R session

Rscript R/clean_chamber_data.R 
Rscript R/clean_dt_data.R
Rscript R/clean_lc_data.R
Rscript R/clean_specimen_data.R
