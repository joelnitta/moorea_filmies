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
    # remove sites on Tahiti (Mt. Aorai)
    filter(str_detect(site, "Aorai", negate = TRUE)),
  
  # - growth habit
  tar_file(filmy_habit_file, "data/filmy_growth_habit.csv"),
  filmy_habit = read_csv(filmy_habit_file),
  
  # - vector of filmy fern species
  filmy_species = filmy_habit$species,
  
  # - phylogenetic tree (time tree, lacks BS values at nodes)
  tar_file(filmy_phy_file, "data/nitta_2017/treepl_Moorea_Tahiti.tre"),
  filmy_phy_no_bs = ape::read.tree(filmy_phy_file) %>%
    # trim to only filmy ferns
    ape::keep.tip(filmy_species),
  
  # - phylogenetic tree (ML tree with BS values at nodes)
  tar_file(filmy_phy_bs_file, "data/nitta_2017/RAxML_bipartitions.all_broad.reduced"),
  filmy_phy_bs = ape::read.tree(filmy_phy_bs_file) %>%
    # trim to only filmy ferns
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
  # remove outliers, sporos measured in field
  light_data = filter_light_data(light_data_all), 
  
  # - specimen data
  tar_file(specimens_raw_file, "data/fern_specimens.csv"),
  specimens_raw = read_csv(specimens_raw_file),
  
  # - gametophyte DT times entered manually (2012 only)
  tar_file(gameto_times_2012_file, "data/2012_gameto_dt_times.csv"),
  gameto_time_summary_2012 = load_gameto_time_summary_2012(gameto_times_2012_file),
  
  # - sporophyte DT times entered manually (2012 only)
  tar_file(sporo_dt_times_file, "data/2012_sporo_dt_times.csv"),
  sporo_dt_times = load_sporo_dt_times(sporo_dt_times_file),
  
  # Process data ----
  
  # - subset collection data to just filmy ferns on Moorea
  filmy_specimens = tidy_filmy_specimens(
    specimens_raw = specimens_raw,
    filmy_species = filmy_species),
  
  # - calculate recovery during DT test per individual
  recovery_indiv = filmy_dt %>%
    # Only include samples that should be presented in the main MS
    filter(section == "main") %>%
    # Filter out individuals with low pre-treatment yields
    filter(!is.na(yield_pre), yield_pre > .4) %>%
    calculate_recovery(),
  
  # - calculate water content during DT test per individual
  rel_water_indiv = calculate_indiv_water(filmy_dt),
  
  # - calculate mean maximum daily VPD per datalogger
  mean_vpd = calculate_mean_max_vpd(climate),
  
  # - calculate mean VPD per species for gametophytes
  mean_vpd_gameto = calculate_mean_vpd_gameto(
    mean_vpd = mean_vpd,
    filmy_specimens = filmy_specimens,
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
    filmy_specimens = filmy_specimens,
    phy = filmy_phy,
    moorea_sites = moorea_sites
  ),
  
  # - combine environmental and range data
  env_range_data = purrr::reduce(
    list(
      mean_vpd_gameto,
      mean_vpd_sporo,
      range_comparison),
    full_join,
    by = "species"
  )  %>% assert(is_uniq, species),
  
  # - fit light curve models to data
  light_models = fit_lc_model(light_data),
  
  # - extract fitted data points
  filmy_lc_model_fitted_data = extract_fitted_lc_data(light_models),
  
  # - extract model parameters: 
  #  ETR at 95% of estimated max value and PAR at that value
  filmy_lc_model_params = extract_lc_model_params(light_models),
  
  # - calculate mean light curve parameters by species and generation
  light_species_means = calculate_mean_light(filmy_lc_model_params),

  # - calculate mean DT recovery by species and generation
  recovery_species_means = calculate_mean_recovery(recovery_indiv),
  
  # - calculate relative water content by species (sporophytes only)
  rel_water_species_means = calculate_mean_water(rel_water_indiv),
  
  # - combine the species means into a single dataframe
  combined_species_means = combine_mean_phys_traits(
    recovery_species_means = recovery_species_means,
    light_species_means = light_species_means
  ),
  
  # - make table of gametophyte DT treatment batches (groups)
  gameto_indiv_times = prepare_gameto_dt_groups(filmy_dt, gameto_time_summary_2012, filmy_specimens),
  
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
  
  # Run PGLS for VPD and range of gameto beyond sporo vs. desiccation tolerance
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
  
  # Track docx style
  tar_file(word_style, "ms/journal-of-plant-research.docx"),
  
  # Render data README
  tar_render(
    data_readme,
    "ms/data_readme.Rmd",
    output_dir = here::here("results/data_readme")
  ),

  # Render MS:
  # - pdf for preprint
  tar_render(
    preprint_pdf, 
    "ms/manuscript.Rmd", 
    output_dir = here::here("results/ms"),
    output_file = "moorea_filmies_preprint.pdf",
    params = list(output_type = "preprint")),
  
  # - pdf for conversion to docx
  tar_render(
    manuscript_pdf, 
    "ms/manuscript.Rmd", 
    output_dir = here::here("results/ms"),
    params = list(output_type = "ms")
    ),
  
  # - docx for submission
  # FIXME: some fixes need to be made manually for MS submission
  # - remove extra space from bib entries
  # - add missing table header rows
  ms_docx = latex2docx(
    latex = "results/ms/manuscript.tex",
    docx = "results/ms/manuscript.docx",
    template = word_style,
    lua_filter = "ms/pagebreak.lua",
    wd = here::here("results"),
    depends = manuscript_pdf
  ),
  
  # Render SI: figures (ESM 1) for submission
  tar_render(
    si_pdf, 
    "ms/si.Rmd", 
    output_file = "ESM_1.pdf", 
    output_dir = here::here("results/si")
    ),
  
  # Render SI: figures (ESM 1) for pre-print
  tar_render(
    si_pdf_bioarxiv, 
    "ms/si.Rmd", 
    output_file = "ESM_1_preprint.pdf", 
    params = list(output_type = "preprint"),
    output_dir = here::here("results/si")
  ),
  
  # Render SI: tables (ESM 2) for sumbission
  tar_file(
    gameto_indiv_times_table_ms,
    write_gameto_times(gameto_indiv_times, ms_type = "ms", "results/si/ESM_2.csv")),
  
  # Render SI: tables (ESM 2) for pre-print
  tar_file(
    gameto_indiv_times_table_preprint,
    write_gameto_times(gameto_indiv_times, ms_type = "preprint", "results/si/ESM_2_preprint.csv"))
  
)