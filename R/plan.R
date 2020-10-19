# Plan the analysis

plan <- drake_plan(

  # Load raw data ----
  
  # Unzip data from Nitta et al. 2017 Ecol. Monographs from Dryad
  # The dataset must be downloaded first by going to
  # https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g, clicking on
  # "Download dataset", and saving to the "data" folder in this project.
  nitta_2017_data = unzip_nitta_2017(
    zipped_path = file_in("data/doi_10.5061_dryad.df59g__v1.zip"),
    unzip_path = "data/nitta_2017",
    # Track data files used as input in analyses
    out1 = file_out("data/nitta_2017/sites.csv"),
    out2 = file_out("data/nitta_2017/treepl_Moorea_Tahiti.tre"),
    out3 = file_out("data/nitta_2017/all_plots.csv")),

  # -microclimate data
  # hobo = readr::read_csv("data/filmy_hobo_data.csv"),

  # - site data
  moorea_sites = readr::read_csv(
    file_in("data/nitta_2017/sites.csv"),
    col_types = "cnnn") %>%
    filter(str_detect(site, "Aorai", negate = TRUE)),

  # -traits
  traits = readr::read_csv("data/filmy_trait_data.csv") %>%
    rename(species = genus_species),

  filmy_species = traits$species,
  
  # -phylogenetic tree
  filmy_phy = ape::read.tree(
    file_in("data/nitta_2017/treepl_Moorea_Tahiti.tre"
    )) %>%
    ape::keep.tip(filmy_species),

  # -physiological data
  recovery_data_raw = readr::read_csv("data/filmy_DT_data.csv"),

  # - community data
  community_matrix_raw = readr::read_csv("data/nitta_2017/all_plots.csv"),

  # - determine which species have widespread gametophytes
  range_types = analyze_dist_pattern(
    community_matrix_raw = community_matrix_raw,
    filmy_species = filmy_species,
    phy = filmy_phy,
    moorea_sites = moorea_sites
  ),

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
  combined_species_means = combine_mean_phys_traits(
    recovery_species_means = recovery_species_means,
    etr_species_means = etr_species_means,
    par_species_means = par_species_means
  ),

  # Phylogenetic signal ----
  phylosig = analyze_phylosig_by_generation(
    combined_species_means = combined_species_means,
    phy = filmy_phy,
    traits_select = c("recover_mean", "etr_mean", "par_mean")
  ),
  
  # Generalized linear mixed models ----
  
  # Uses phylogeny for DT only
  glmms = run_glmm(
    combined_species_means = combined_species_means, 
    traits = traits,
    phy = filmy_phy)

)
