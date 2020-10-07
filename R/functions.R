# Define functions used in the analysis

#' Unzip Nitta et al 2017 Ecol Mono data file downloaded
#' from Dryad (https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g)
#' and extract needed data files.
#' 
#' (This is used to download the file as part of the function, but
#' as it is not clear how to do this with the new Dryad API that
#' was implemented ca. Sept 2019, it now requires the user to have 
#' already downloaded the zipped data file from Dryad manually.)
#' 
#' @param zipped_path Name of downloaded zip file.
#' @param unzip_path Path to directory to put the unzipped
#' contents (will be created if needed).
#' @param ... Extra arguments; not used by this function, but
#' meant for tracking with drake.
#' @return Four unzipped data files:
#' - all_plots.csv: Community matrix for ferns of Moorea and Tahiti
#' including sporophytes (plot names with "_S") and gametophytes
#' (plot names with "_G").
#' - treepl_Moorea_Tahiti.tre: Dated tree for pteridophytes of Moorea and Tahiti
#' - sites.csv: Site metadata including name, latitude, longitude, and elevation (m)
#' - species.csv: List of species included in the Nitta et al. 2017 study
#'
unzip_nitta_2017 <- function (zipped_path, unzip_path, ...) {
  
  # The dryad file is actually a zipped folder within a zipped folder.
  # Unzip the first one to a temporary folder, then unzip the rest from there.
  temp_dir <- tempdir()
  
  unzip(zipped_path, exdir = temp_dir, overwrite = TRUE)
  
  temp_zip <- fs::path(temp_dir, "data_and_scripts.zip")
  
  # Unzip only needed data files to data/nitta_2017/
  unzip(temp_zip, "data_and_scripts/Comm_Phylo_Analysis/data/all_plots.csv", exdir = unzip_path, junkpaths = TRUE, overwrite = TRUE)
  unzip(temp_zip, "data_and_scripts/Comm_Phylo_Analysis/data/treepl_Moorea_Tahiti.tre", exdir = unzip_path, junkpaths = TRUE, overwrite = TRUE)
  unzip(temp_zip, "data_and_scripts/shared_data/sites.csv", exdir = unzip_path, junkpaths = TRUE, overwrite = TRUE)
  
  # Cleanup
  fs::file_delete(temp_zip)
  
}

#' Process raw community matrix 
#'
#' @param community_matrix_path Path to community matrix of Nitta et al 2017
#' @param species_list Vector of species to include
#' @param moorea_sites Dataframe of sites on Moorea in this study
#'
#' @return Tibble (community data matrix)
#' 
process_community_matrix <- function (community_matrix_path, species_list, moorea_sites) {
  read_csv(community_matrix_path) %>%
    rename(site = X1) %>% 
    gather(species, abundance, -site) %>%
    filter(str_detect(site, "_S")) %>%
    mutate(site = str_remove_all(site, "_S")) %>%
    # Keep only sites on Moorea
    filter(site %in% moorea_sites$site) %>%
    # Keep only species in species list (ferns of Moorea)
    filter(species %in% species_list) %>%
    spread(site, abundance)
}

clean_recovery <- function (data) {
  
  # Convert any recovery >100% to 1
  data %>%
    mutate(recover = case_when(
      recover > 1 ~ 1,
      TRUE ~ recover
    )) %>%
    rename(
      specific_epithet = species,
      species = genus_species)
  
}

clean_light <- function (data) {
  
  # Exclude outliers
  data %>%
    filter(exclude_sample == 0) %>%
    filter(exclude_points == 0) %>%
    rename(species = genus_species)
  
}

calculate_mean_recovery <- function (data) {
  
  # Calculate mean recovery by salt treatment, dry time, recovery time, generation, and species
  data %>%
    group_by(salt, drytime, rectime, generation, species) %>%
    summarize(
      recover_mean = mean(recover, na.rm = TRUE),
      recover_sd = sd(recover, na.rm = TRUE),
      recover_n = n()
    ) %>%
    ungroup()
  
}

calculate_mean_etr <- function (data) {
  
  individual_etr <-
    data %>%
    assert(not_na, ETR) %>%
    group_by(species, generation, individual) %>%
    summarize(
      etr_max = max(ETR)
    ) %>%
    ungroup
  
  # Next calculate the mean value by species and generation
  individual_etr %>%
    group_by(generation, species) %>%
    summarize(
      etr_mean = mean(etr_max),
      etr_sd = sd(etr_max),
      etr_n = n()
    ) %>%
    ungroup()
  
}

calculate_mean_par <- function (data) {
  
  # First calculate critical PAR value for each individual
  individual_par <-
    data %>%
    select(species, generation, individual, PAR, ETR) %>%
    nest(data = c(PAR, ETR)) %>%
    mutate(
      nls_mod = map(
        data, 
        # nonlinear least-squares estimates of parameters of the subsetted data
        ~nls(
          ETR~max(ETR)*(1-exp(-k*PAR)),
          start = list(k = 0.04),
          data =.,
          trace = FALSE,
          control = list(maxiter = 500))
      ),
      k_stats = map(nls_mod, broom::tidy)
    ) %>%
    select(species, generation, individual, k_stats) %>%
    unnest(cols = c(k_stats)) %>%
    # calculate critical PAR (PAR where reach 95% of max ETR)
    mutate(par_critical = -log(0.05)/estimate)
  
  # Next calculate the mean value by species and generation
  individual_par %>%
    group_by(generation, species) %>%
    summarize(
      par_mean = mean(par_critical),
      par_sd = sd(par_critical),
      par_n = n()
    ) %>%
    ungroup()
  
}

combine_mean_phys_traits <- function (recovery_species_means, etr_species_means, par_species_means) {
  
  recover_filtered <-
    recovery_species_means %>%
    mutate(
      keep = case_when(
        species == "Callistopteris_apiifolia" & generation == "sporo" & salt == "NaCl" & drytime == "2d" & rectime == "48hr" ~ TRUE,
        species == "Abrodictyum_dentatum" & generation == "sporo" & salt == "NaCl" & drytime == "2d" & rectime == "48hr" ~ TRUE,
        salt == "MgNO3" & drytime == "2d" & rectime == "48hr" ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(keep) %>%
    select(species, generation, starts_with("recover"))
  
  list(
    recover_filtered, etr_species_means, par_species_means
  ) %>%
    reduce(left_join, by = c("species", "generation"))
  
}


# Phylogentic signal ----

#' Analyze phylogenetic signal in a continuous trait of interest
#'
#' @param selected_trait Name of trait to analyze phylogenetic signal
#' @param traits Dataframe including all untransformed traits, with
#' 'species' as a column and other traits as other columns.
#' @param phy Phylogeny
#'
#' @return List of estimated Blomberg's K and Pagel's lambda and
#' their significance
#' 
analyze_cont_phylosig <- function (selected_trait, traits, phy) {
  
  # Trim data to non-missing trait values and
  # make sure species in same order in tree and traits
  traits_select <- traits %>% select(species, all_of(selected_trait)) %>%
    remove_missing(na.rm = TRUE)
  
  traits_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "traits") 
  phy_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "tree") 
  
  # Extract named vector of trait values for phylosig()
  trait_vec <- pull(traits_trim, selected_trait) %>%
    set_names(traits_trim$species)
  
  # Run phylosig() on selected trait
  # using Blomberg's K
  k_model <- phytools::phylosig(phy_trim, trait_vec, method = "K", test = TRUE)
  # and Pagel's lambda
  lambda_model <- phytools::phylosig(phy_trim, trait_vec, method = "lambda", test = TRUE)
  
  # get model results
  list(trait = selected_trait,
       kval = k_model$K,
       k_pval = k_model$P,
       lambda = lambda_model$lambda,
       lambda_pval = lambda_model$P)
  
}

#' Analyze phylogenetic signal by generation (sporophyte or gametophyte)
#'
#' @param combined_species_means Combined species means, including a 'species' column, columns for each trait,
#' and 'generation' column (either 'sporo' or 'gameto')
#' @param phy Phylogeny
#' @param traits_select Character vector of traits to analyze (must match column names in combined_species_means)
#'
#' @return Dataframe
#' 
analyze_phylosig_by_generation <- function(combined_species_means, phy, traits_select = c("recover_mean", "etr_mean", "par_mean")) {
  
  # Subset to only gametophytes
  gameto_traits <- combined_species_means %>%
    filter(generation == "gameto")
  
  # Analyze phylogenetic signal
  gameto_phylosig <- map_df(traits_select, ~analyze_cont_phylosig(., gameto_traits, phy)) %>%
    mutate(generation = "gameto")
  
  # Subset to only sporophytes
  sporo_traits <- combined_species_means %>%
    filter(generation == "sporo")
  
  # Analyze phylogenetic signal
  sporo_phylosig <- map_df(traits_select, ~analyze_cont_phylosig(., sporo_traits, phy)) %>%
    mutate(generation = "sporo")
  
  # Combine results
  bind_rows(sporo_phylosig, gameto_phylosig)
  
}

# Etc ----

#' Match trait data and tree
#' 
#' Order of species in traits will be rearranged to match the
#' phylogeny.
#'
#' @param traits Dataframe of traits, with 'species' column and
#' additional columns, one for each trait
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the traits, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
match_traits_and_tree <- function (traits, phy, return = c("traits", "tree")) {
  
  assert_that("species" %in% colnames(traits))
  
  # Keep only species in phylogeny
  traits <- traits %>%
    filter(species %in% phy$tip.label) 
  
  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, traits$species))
  
  # Get traits in same order as tips
  traits <- left_join(
    tibble(species = phy$tip.label),
    traits,
    by = "species"
  )
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(traits$species, phy$tip.label)))
  
  # Return traits or tree
  assert_that(return %in% c("tree", "traits"))
  
  if(return == "tree") { 
    return (phy) 
  } else {
    return (traits)
  }
  
}

#' Match community data and tree
#' 
#' Order of species in comm will be rearranged to match the
#' phylogeny.
#'
#' @param comm Community data frame, with one column for sites and
#' the rest for species.
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the community, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
match_comm_and_tree <- function (comm, phy, return = c("comm", "tree")) {
  
  assert_that("species" %in% colnames(comm))
  
  # Keep only species in phylogeny
  comm <- comm %>%
    filter(species %in% phy$tip.label) 
  
  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, comm$species))
  
  # Get comm in same order as tips
  comm <- left_join(
    tibble(species = phy$tip.label),
    comm
  )
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, phy$tip.label)))
  
  # Return comm or tree
  assert_that(return %in% c("tree", "comm"))
  
  if(return == "tree") { 
    return (phy) 
  } else {
    return (comm)
  }
  
}

# Plots ----

make_sporo_dt_plot <- function  (recovery_species_means, traits) {
  # Define color-blind palette
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  plot_data <-
    recovery_species_means %>%
    # Use only sporophytes
    filter(generation == "sporo") %>%
    # Don't plot control
    filter(salt != "Control") %>%
    # Add growth habit data
    left_join(traits, by = "species") %>%
    assert(not_na, habit, range) %>%
    # Abbreviate species names to G. species
    mutate(species = str_replace_all(species, "[a-z]+_", ". ")) %>%
    # Reorder levels for species: alphabetical within growth habit
    mutate(species = fct_reorder2(species, habit, species, .desc = FALSE)) %>%
    rename(recover = recover_mean, sd = recover_sd, n = recover_n)
  
  # Self-defined formatting function for percent
  percent_formatter <- function(x) {
    lab <- x*100
  }
  
  # Format legend (include title, but get rid of boxes)
  legend_theme <- theme (legend.background = element_rect(colour = "transparent", fill = "transparent", size = 0.5),
                         legend.key = element_rect(colour = "white", fill = "white", size = 0.5),
                         legend.text = element_text(size = 12))
  
  # Set colors for desiccation intensity
  col <- c(cbPalette[7], cbPalette[5], cbPalette[6])
  
  # Set offset of datapoints (dodge)
  pd <- position_dodge(.3)
  
  # Make plot
  ggplot(data = plot_data, aes(x = rectime, y = recover, group = interaction(drytime, salt), shape = salt, fill = salt)) +
    geom_line(position = pd, aes(linetype = drytime)) +
    geom_errorbar(
      aes(ymin = recover - sd, ymax = recover + sd), 
      colour = "grey50", width = 0, 
      position = pd) +
    geom_point(position = pd, size = 1.5) +
    scale_linetype_manual(name = "Desiccation Time",
                          values = c("solid", "dashed"),
                          breaks = c("15d", "2d"),
                          labels = c("15 days", "2 days")) +
    scale_shape_manual(name = "Desiccation Intensity",
                       values = c(23,24,21),
                       breaks = c("NaCl", "MgNO3",  "LiCl"),
                       labels = c("-38 MPa", "-86 MPa", "-282 Mpa")) + 
    scale_fill_manual(name = "Desiccation Intensity",
                      values = col, 
                      breaks = c("NaCl", "MgNO3",  "LiCl"),
                      labels = c("-38 MPa", "-86 MPa", "-282 Mpa")) + 
    scale_y_continuous(label = percent_formatter, limits = c(0,1.1), expand = c(0,0), breaks=c(0,.2,.4,.6,.8,1.0)) +
    scale_x_discrete(limits = c("30min","24hr","48hr"),
                     labels = c("0.5", "24", "48")) +
    xlab("Recovery Time (hr)") +
    ylab("Recovery (%)") +
    theme_bw() +
    legend_theme +
    theme(legend.position="bottom",
          panel.grid.minor=element_blank(), 
          panel.grid.major=element_blank()) + 
    facet_wrap( ~ species, ncol=3)
}

#' Determine which species have widespread gametophytes
#'
#' @param raw_community_matrix Community matrix of Nitta et al 2017
#' @param filmy_species Character vector of filmy fern species names
#' @param phy Phylogeny including filmy ferns of Moorea and Tahiti
#' @param moorea_sites Site data for plots on Moorea including elevation
#' 
#' @return Dataframe
#' 
analyze_dist_pattern <- function(community_matrix_raw, filmy_species, phy, moorea_sites) {
  
  # Convert raw community data to long format by sporophyte and gametophyte
  # min and max elevational range.
  range_long <-
    community_matrix_raw %>%
    rename(site = X1) %>% 
    gather(species, abundance, -site) %>%
    filter(species %in% filmy_species) %>%
    mutate(generation = case_when(
      str_detect(site, "_S") ~ "sporophyte",
      str_detect(site, "_G") ~ "gametophyte"
    )) %>%
    mutate(site = str_remove_all(site, "_S") %>% str_remove_all("_G")) %>%
    # Keep only sites on Moorea
    filter(site %in% moorea_sites$site) %>%
    left_join(moorea_sites, by = "site") %>%
    filter(abundance > 0) %>%
    assert(not_na, el) %>%
    group_by(species, generation) %>%
    summarize(
      min_range = min(el),
      max_range = max(el)
    ) %>%
    ungroup
  
  # Compare ranges to see which species have gametophytes beyond sporophytes
  full_join(
    range_long %>%
      gather(variable, value, -species, -generation) %>%
      filter(generation == "sporophyte") %>% 
      select(-generation) %>%
      spread(variable, value) %>%
      rename_at(vars(max_range, min_range), ~paste("sporo", ., sep = "_")),
    range_long %>%
      gather(variable, value, -species, -generation) %>%
      filter(generation == "gametophyte") %>% 
      select(-generation) %>%
      spread(variable, value) %>%
      rename_at(vars(max_range, min_range), ~paste("gameto", ., sep = "_")),
    by = "species"
  ) %>%
    assert(is_uniq, species) %>%
    mutate(
      range = case_when(
        gameto_max_range > sporo_max_range + 200 ~ "widespread",
        gameto_min_range < sporo_min_range - 200 ~ "widespread",
        is.na(sporo_max_range) & is.na(sporo_min_range) ~ "widespread",
        TRUE ~ "not_widespread"
      )
    )
  
}
