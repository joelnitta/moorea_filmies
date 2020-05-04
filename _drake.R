# Load packages, functions, and plan
source("R/packages.R")
source("R/functions.R")
source("R/plan.R")

# Set up drake configuration
drake_config(plan, verbose = 1, seed = 0)
