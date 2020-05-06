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
    ))
  
}

clean_light <- function (data) {
  
  # Exclude outliers
  data %>%
    filter(exclude_sample == 0) %>%
    filter(exclude_points == 0)
  
}

calculate_mean_recovery <- function (data) {
  
  # Calculate mean recovery by salt treatment, dry time, recovery time, generation, and species
  data %>%
    group_by(salt, drytime, rectime, generation, genus_species) %>%
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
    group_by(genus_species, generation, individual) %>%
    summarize(
      etr_max = max(ETR)
    ) %>%
    ungroup
  
  # Next calculate the mean value by species and generation
  individual_etr %>%
    group_by(generation, genus_species) %>%
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
    select(genus_species, generation, individual, PAR, ETR) %>%
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
    select(genus_species, generation, individual, k_stats) %>%
    unnest(cols = c(k_stats)) %>%
    # calculate critical PAR (PAR where reach 95% of max ETR)
    mutate(par_critical = -log(0.05)/estimate)
  
  # Next calculate the mean value by species and generation
  individual_par %>%
    group_by(generation, genus_species) %>%
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
        genus_species == "Callistopteris_apiifolia" & generation == "sporo" & salt == "NaCl" & drytime == "2d" & rectime == "48hr" ~ TRUE,
        genus_species == "Abrodictyum_dentatum" & generation == "sporo" & salt == "NaCl" & drytime == "2d" & rectime == "48hr" ~ TRUE,
        salt == "MgNO3" & drytime == "2d" & rectime == "48hr" ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(keep) %>%
    select(genus_species, generation, starts_with("recover"))
  
  list(
    recover_filtered, etr_species_means, par_species_means
  ) %>%
    reduce(left_join, by = c("genus_species", "generation"))
  
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
