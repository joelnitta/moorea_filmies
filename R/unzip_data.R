source("R/packages.R")
source("R/functions.R")

# Unzip data from Nitta et al. 2017 Ecol. Monographs from Dryad
# The dataset must be downloaded first by going to
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g, clicking on
# "Download dataset", and saving to the "data" folder in this project.
unzip_nitta_2017(
  zipped_path = "data/doi_10.5061_dryad.df59g__v1.zip",
  unzip_path = "data/nitta_2017")

# Unzip data from Nitta et al. 2020 New Phyt. from Dryad
# The dataset must be downloaded first by going to
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.fqz612jps, clicking on
# "Download dataset", and saving to the "data" folder in this project.
unzip(
  zipfile = "data/doi_10.5061_dryad.fqz612jps__v2.zip", 
  files = "moorea_climate.csv", 
  exdir = "data/nitta_2020/", 
  junkpaths = TRUE, 
  overwrite = TRUE)

# Unzip data from this paper
# The dataset must be downloaded first by going to
# https://figshare.com/s/d6349abf01a3756a5aae, clicking on
# "Download all", and saving to the "data" folder in this project.
unzip(
  zipfile = "data/14184572.zip", 
  exdir = "data/", 
  junkpaths = TRUE, 
  overwrite = TRUE)
