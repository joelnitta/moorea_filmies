#!/usr/bin/bash

# Run each data cleaning script in its own R session

# Data files "filmy_growth_habit.csv" and those extracted by R/unzip_data.R must
# already be present in data/

Rscript R/clean_specimen_data.R
Rscript R/clean_chamber_data.R 
Rscript R/clean_dt_data.R
Rscript R/clean_lc_data.R
