source("R/packages.R")
source("R/functions.R")

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
  
  filmy_species = filmy_habit$species,
  
  # - phylogenetic tree
  tar_file(filmy_phy_file, "data/nitta_2017/treepl_Moorea_Tahiti.tre"),
  filmy_phy = ape::read.tree(filmy_phy_file) %>%
    ape::keep.tip(filmy_species),
  
  # - sporophyte desiccation tolerance
  tar_file(filmy_sporo_dt_file, "data/filmy_sporo_dt.csv"),
  filmy_sporo_dt = load_sporo_dt(filmy_sporo_dt_file),
  
  # - gametophyte desiccation tolerance
  tar_file(filmy_gameto_dt_file, "data/filmy_gameto_dt.csv"),
  filmy_gameto_dt = load_gameto_dt(filmy_gameto_dt_file),
  
  # - community data
  tar_file(community_matrix_raw_file, "data/nitta_2017/all_plots.csv"),
  community_matrix_raw = read_csv(community_matrix_raw_file),

  # - light responses
  tar_file(light_data_file, "data/filmy_light_curves.csv"),
  light_data = read_csv(light_data_file),
  
  # - specimen data
  tar_file(specimens_raw_file, "data/fern_specimens.csv"),
  specimens_raw = read_csv(specimens_raw_file),
  
  # Process data ----
  
  # - calculate recovery during DT test per individual
  gameto_recovery_indiv = filmy_gameto_dt %>%
    # Filter out individuals with low pre-treatment yields
    filter(yield_pre > 400) %>%
    calc_recovery(),
  
  sporo_recovery_indiv = filmy_sporo_dt %>%
    # Only include samples that should be presented in the main MS
    filter(section == "main") %>%
    # Filter out individuals with low pre-treatment yields
    filter(yield_pre > 400) %>%
    calc_recovery(),
  
  recovery_indiv = bind_rows(gameto_recovery_indiv, sporo_recovery_indiv),
  
  # - calculate mean VPD per datalogger
  mean_vpd = calculate_mean_vpd(climate),
  
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
  
  # - determine which species have widespread gametophytes
  range_types = analyze_dist_pattern(
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
      range_types),
    left_join,
    by = "species"
  ),
  
  # - calculate PAR95 by individual
  par_indiv = calculate_indiv_par(light_data),
  
  # - calculate ETRmax by individual
  etr_indiv = calculate_indiv_etr(light_data),
  
  # - calculate mean DT recovery by species and generation
  recovery_species_means = calc_mean_recovery(recovery_indiv),
  
  # - calculate mean ETR max by species and generation
  etr_species_means = calculate_mean_etr(etr_indiv),
  
  # - calculate mean PAR95 by species and generation
  par_species_means = calculate_mean_par(par_indiv),
  
  # - combine the species means into a single df
  combined_species_means = combine_mean_phys_traits(
    recovery_species_means = recovery_species_means,
    etr_species_means = etr_species_means,
    par_species_means = par_species_means
  ),
  
  # t-test ----
  
  # Perform two-sided t-test on DT and light responses across generations
  t_test_results = run_t_test(
    recovery_data = recovery_indiv, 
    par_indiv = par_indiv, 
    etr_indiv = etr_indiv),
  
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
  
  # PGLS ----
  # (Phylogenetic Generalized Least Squares)
  
  # Run PGLS for gametophyte range size vs. desiccation tolerance
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
  tar_render(manuscript_pdf, "ms/manuscript.Rmd", output_dir = here::here("results")),
  
  # Render MS Word
  ms_docx = latex2docx(
    latex = "results/manuscript.tex",
    docx = "results/manuscript.docx",
    template = "ms/journal-of-plant-research.docx",
    wd = here::here("results"),
    depends = manuscript_pdf
  )
  
)