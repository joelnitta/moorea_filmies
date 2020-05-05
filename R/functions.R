# Define functions used in the analysis

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
