# Define functions used in the analysis

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

