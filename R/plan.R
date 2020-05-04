# Plan the analysis

plan <- drake_plan(
  
  # Load raw data ----
  
  # -microclimate data
  hobo = readr::read_csv("data/filmy_hobo_data.csv"),
  
  # -site data
  lat_long_el = readr::read_csv("data/filmy_sites.csv"),
  
  # -phylogenetic tree
  phy = ape::read.tree("data/filmy_tree.tre"),
  
  # -traits
  traits = readr::read_csv("data/filmy_trait_data.csv"),
  
  # -physiological data
  recovery_data_raw = readr::read_csv("data/filmy_DT_data.csv"),

  # FIXME: check number of measurements per individual; shouldn't have more than 10 each.
  # Crepidomanes_minutum2 (2834), Crepidomanes_minutum2 (2998), and Hymenophyllum_braithwaitei (Hymenophyllum_sp1_6)
  # have too many
  light_data_raw = readr::read_csv("data/filmy_light_data.csv"),
  
  # Clean raw data ----
  recovery_data = clean_recovery(recovery_data_raw),
  
  light_data = clean_light(light_data_raw),
  
  # Calculate means by species and generation ----
  # -% recovery
  recovery_species_means = calculate_mean_recovery(recovery_data),
  
  # -ETR max
  etr_species_means = calculate_mean_etr(light_data),
  
  # -PAR 95%
  par_species_means = calculate_mean_par(light_data),
  
  # -Combine the means into a single df
  combined_species_means = combine_mean_phys_traits
  (recovery_species_means = recovery_species_means, 
    etr_species_means = etr_species_means, 
    par_species_means = par_species_means),
  
)
