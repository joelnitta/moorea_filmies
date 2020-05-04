# Plan the analysis

plan <- drake_plan(
  
  # Load data ----
  
  # Raw microclimate data
  hobo = readr::read_csv("data/filmy_hobo_data.csv"),
  
  # Site data
  lat_long_el = readr::read_csv("data/filmy_sites.csv"),
  
  # Phylogenetic tree
  phy = ape::read.tree("data/filmy_tree.tre"),
  
  # Traits
  traits = readr::read_csv("data/filmy_trait_data.csv"),
  
  # Physiological data
  DT_data = readr::read_csv("data/filmy_DT_data.csv"),
  light_data = readr::read_csv("data/filmy_light_data.csv"),
  
  # Presence/absence community data
  comm.s = readr::read_csv("data/filmy_sporo_comm.csv"),
  comm.g = readr::read_csv("data/filmy_gameto_comm.csv")
  
  
)