source("R/packages.R")
source("R/functions.R")

# Define workflow plan
tar_plan(
  
  # Load data ----
  
  # - microclimate data
  tar_file(climate_file, "data/nitta_2020/moorea_climate.csv"),
  climate = read_csv(climate_file),
  
  # - site data
  tar_file(moorea_sites_file, "data/nitta_2017/sites.csv"),
  moorea_sites = read_csv(
    moorea_sites_file,
    col_types = "cnnn") %>%
    filter(str_detect(site, "Aorai", negate = TRUE)),
  
  # - growth habit
  tar_file(filmy_habit_file, "data/filmy_growth_habit.csv"),
  filmy_habit = read_csv(filmy_habit_file),
  
  # - vector of filmy fern species
  filmy_species = filmy_habit$species,
  
  # - phylogenetic tree (time tree, lacks BS values at nodes)
  tar_file(filmy_phy_file, "data/nitta_2017/treepl_Moorea_Tahiti.tre"),
  filmy_phy_no_bs = ape::read.tree(filmy_phy_file) %>%
    ape::keep.tip(filmy_species),
  
  # - phylogenetic tree (ML tree with BS values at nodes)
  tar_file(filmy_phy_bs_file, "data/nitta_2017/RAxML_bipartitions.all_broad.reduced"),
  filmy_phy_bs = ape::read.tree(filmy_phy_bs_file) %>%
    ape::keep.tip(filmy_phy_no_bs$tip.label),
  
  # - transfer bootstrap values from ML tree to time tree
  filmy_phy = transfer_bs(filmy_phy_bs, filmy_phy_no_bs),
  
  # - desiccation tolerance yields
  tar_file(filmy_dt_file, "data/filmy_dt.csv"),
  filmy_dt = load_filmy_dt(filmy_dt_file),
  
  # - desiccation tolerance chamber temp and RH
  tar_file(filmy_dt_chamber_file, "data/filmy_dt_chamber.csv"),
  filmy_dt_chamber = load_filmy_dt_chamber(filmy_dt_chamber_file),
  
  # - community data
  tar_file(community_matrix_raw_file, "data/nitta_2017/all_plots.csv"),
  community_matrix_raw = read_csv(community_matrix_raw_file),

  # - light responses
  tar_file(light_data_file, "data/filmy_light_curves.csv"),
  light_data_all = load_filmy_lc(light_data_file),
  light_data = filter(light_data_all, outlier == FALSE), # remove outliers
  
  # - specimen data
  tar_file(specimens_raw_file, "data/fern_specimens.csv"),
  specimens_raw = read_csv(specimens_raw_file),
  
  # Process data ----
  
  # - calculate recovery during DT test per individual
  recovery_indiv = filmy_dt %>%
    # Only include samples that should be presented in the main MS
    filter(section == "main") %>%
    # Filter out individuals with low pre-treatment yields
    filter(!is.na(yield_pre), yield_pre > 400) %>%
    calculate_recovery(),
  
  # - calculate water content during DT test per individual
  rel_water_indiv = calculate_indiv_water(filmy_dt),
  
  # - calculate mean maximum daily VPD per datalogger
  mean_vpd = calculate_mean_max_vpd(climate),
  
  # - calculate mean VPD per species for gametophytes
  mean_vpd_gameto = calculate_mean_vpd_gameto(
    mean_vpd = mean_vpd,
    specimens_raw = specimens_raw,
    moorea_sites = moorea_sites,
    filmy_species = filmy_species
  ),
  
  # - calculate mean VPD per species for sporophytes
  mean_vpd_sporo = calculate_mean_vpd_sporo(
    community_matrix_raw = community_matrix_raw,
    mean_vpd = mean_vpd,
    traits = filmy_habit,
    moorea_sites = moorea_sites,
    filmy_species = filmy_species
  ),
  
  # - determine range of gametophytes beyond sporophytes
  range_comparison = analyze_dist_pattern(
    community_matrix_raw = community_matrix_raw,
    filmy_species = filmy_species,
    phy = filmy_phy,
    moorea_sites = moorea_sites
  ),
  
  # - combine environmental and range data
  env_range_data = purrr::reduce(
    list(
      mean_vpd_gameto,
      mean_vpd_sporo,
      range_comparison),
    left_join,
    by = "species"
  ),
  
  # - fit light curve models to data
  light_models = fit_lc_model(light_data),
  
  # - extract model parameters: 
  # critical PAR, and ETR at 95% of estimated max value
  filmy_lc_model_params = extract_lc_model_params(light_models),
  
  # - calculate mean light curve parameters by species and generation
  light_species_means = calculate_mean_light(filmy_lc_model_params),

  # - calculate mean DT recovery by species and generation
  recovery_species_means = calculate_mean_recovery(recovery_indiv),
  
  # - calculate relative water content by species (sporophytes only)
  rel_water_species_means = calculate_mean_water(rel_water_indiv),
  
  # - combine the species means into a single df
  combined_species_means = combine_mean_phys_traits(
    recovery_species_means = recovery_species_means,
    light_species_means = light_species_means
  ),
  
  # t-test ----
  # Perform two-sided t-test on DT and light responses across generations
  t_test_results = run_t_test(
    filmy_lc_model_params = filmy_lc_model_params, 
    recovery_indiv = recovery_indiv),
  
  # Phylogenetic signal ----
  phylosig = analyze_phylosig_by_generation(
    combined_species_means = combined_species_means,
    phy = filmy_phy,
    traits_select = c("recovery_mean", "etr_mean", "par_mean")
  ),
  
  # GLMMS ----
  # (Generalized Linear Mixed Models)
  
  # Uses phylogeny for DT only
  glmms = run_glmm(
    combined_species_means = combined_species_means, 
    traits = filmy_habit,
    phy = filmy_phy),
  
  # Get DIC for each GLMM
  glmm_summary = tidy_glmms(glmms),
  
  # Get parameters for best GLMMs only
  glmm_params = tidy_best_glmm_params(glmms),
  
  # PGLS ----
  # (Phylogenetic Generalized Least Squares)
  
  # Run PGLS for VPD and gametophyte range breadth vs. desiccation tolerance
  env_range_recover_data = combine_env_env_range_recover(
    combined_species_means = combined_species_means,
    env_range_data = env_range_data
  ),
  
  env_range_dt_model = run_pgls(
    env_range_recover_data = env_range_recover_data,
    phy = filmy_phy),
  
  env_range_dt_model_summary = map_df(env_range_dt_model, tidy_pgls, .id = "model"),
  
  # Render manuscript ----
  
  # Track bibliography files
  tar_file(refs, "ms/references.bib"),
  tar_file(refs_other, "ms/references_other.yaml"),

  # Render PDF
  tar_render(manuscript_pdf, "ms/manuscript.Rmd", output_dir = here::here("results/ms")),
  
  # Render MS Word
  ms_docx = latex2docx(
    latex = "results/ms/manuscript.tex",
    docx = "results/ms/manuscript.docx",
    template = "ms/journal-of-plant-research.docx",
    lua_filter = "ms/pagebreak.lua",
    wd = here::here("results"),
    depends = manuscript_pdf
  ),
  
  # SI
  tar_render(si_pdf, "ms/si.Rmd", output_dir = here::here("results/ms"))
  
)