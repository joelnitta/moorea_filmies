# Plan the analysis

plan <- drake_plan(

  # Unzip raw data ----
  
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

  # Unzip data from Nitta et al. 2020 New Phyt. from Dryad
  # The dataset must be downloaded first by going to
  # https://datadryad.org/stash/dataset/doi:10.5061/dryad.fqz612jps, clicking on
  # "Download dataset", and saving to the "data" folder in this project.
  climate_unzipped = unzip(
    zipfile = file_in("data/doi_10.5061_dryad.fqz612jps__v4.zip"), 
    files = "moorea_climate.csv", 
    exdir = file_out("data/nitta_2020/"), 
    junkpaths = TRUE, 
    overwrite = TRUE),
  
  # Load data ----
  
  # - microclimate data
  climate = readr::read_csv(file_in("data/nitta_2020/moorea_climate.csv")),

  # - site data
  moorea_sites = readr::read_csv(
    file_in("data/nitta_2017/sites.csv"),
    col_types = "cnnn") %>%
    filter(str_detect(site, "Aorai", negate = TRUE)),

  # - growth habit
  filmy_habit = readr::read_csv("data/filmy_growth_habit.csv"),

  filmy_species = filmy_habit$species,
  
  # - phylogenetic tree
  filmy_phy = ape::read.tree(
    file_in("data/nitta_2017/treepl_Moorea_Tahiti.tre"
    )) %>%
    ape::keep.tip(filmy_species),
  
  # - sporophyte desiccation tolerance
  filmy_sporo_dt = load_sporo_dt("data/filmy_sporo_dt.csv"),
  
  # - gametophyte desiccation tolerance
  filmy_gameto_dt = load_gameto_dt("data/filmy_gameto_dt.csv"),
  
  # - physiological data
  recovery_data_raw = readr::read_csv("data/filmy_DT_data.csv"),

  # - community data
  community_matrix_raw = readr::read_csv("data/nitta_2017/all_plots.csv"),

  # FIXME: check number of measurements per individual; shouldn't have more than 10 each.
  # Crepidomanes_minutum2 (2834), Crepidomanes_minutum2 (2998), and Hymenophyllum_braithwaitei (Hymenophyllum_sp1_6)
  # have too many
  light_data_raw = readr::read_csv("data/filmy_light_data.csv"),

  # Raw specimen data from specimens spreadsheet
  specimens_raw = read_csv("data/specimens.csv"),
  
  # Process data ----
  filmy_gameto_recovery = calc_recovery(filmy_gameto_dt),
  
  filmy_sporo_recovery =
    filmy_sporo_dt %>%
    # Only include samples that should be presented in the main MS
    filter(section == "main") %>%
    calc_recovery(),
  
  filmy_sporo_recovery_mean = calc_mean_recovery(filmy_sporo_recovery),
  
  filmy_gameto_recovery_mean = calc_mean_recovery(filmy_gameto_recovery),
  
  # FIXME: delete unused functions
  # recovery_data = clean_recovery(recovery_data_raw),

  # FIXME: delete unused functions
  # light_data = clean_light(light_data_raw),
  
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

  # Means by species and generation ----
  # - DT recovery
  recovery_species_means = calculate_mean_recovery(recovery_data),

  # - ETR max
  etr_species_means = calculate_mean_etr(etr_indiv),

  # - PAR95
  par_species_means = calculate_mean_par(par_indiv),

  # - Combine the means into a single df
  combined_species_means = combine_mean_phys_traits(
    recovery_species_means = recovery_species_means,
    etr_species_means = etr_species_means,
    par_species_means = par_species_means
  ),
  
  # t-test ----
  
  # Perform two-sided t-test on DT and light responses across generations
  t_test_results = run_t_test(
    recovery_data = recovery_data, 
    par_indiv = par_indiv, 
    etr_indiv = etr_indiv),

  # Phylogenetic signal ----
  phylosig = analyze_phylosig_by_generation(
    combined_species_means = combined_species_means,
    phy = filmy_phy,
    traits_select = c("recover_mean", "etr_mean", "par_mean")
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
  
  range_dt_model = run_pgls(
    env_range_recover_data = env_range_recover_data, 
    phy = filmy_phy),
  
  # Render manuscript ----
  
  # Track bibliography files
  refs = target("ms/references.bib", format = "file"),
  refs_other = target("ms/references_other.yaml", format = "file"),
  
  # First render to PDF, keeping the latex
  ms_pdf = render_tracked(
    input = knitr_in("ms/manuscript.Rmd"),
    quiet = TRUE,
    output_dir = here::here("results"),
    tracked_output = file_out("results/manuscript.tex"),
    dep1 = refs,
    dep2 = refs_other
  ),
  
  # Next use the latex to convert to docx with pandoc
  ms_docx = latex2docx(
    latex = file_in("results/manuscript.tex"),
    docx = file_out("results/manuscript.docx"),
    template = file_in("ms/journal-of-plant-research.docx"),
    wd = here::here("results")
  )

)
