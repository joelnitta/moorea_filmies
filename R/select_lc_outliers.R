# Select outliers in light curve data

# Load packages
source(here::here("R/packages.R"))

# Load functions
source("R/functions.R")

# Load raw light-curve data
filmy_lc <- read_csv("data/filmy_light_curves.csv", col_types = cols(
  type = col_character(),
  no = col_double(),
  f = col_double(),
  fm = col_double(),
  par = col_double(),
  yield = col_double(),
  etr = col_double(),
  date_time = col_datetime(format = ""),
  species = col_character(),
  individual = col_character(),
  generation = col_character(),
  coll_num = col_character(),
  sporo_site = col_character(),
  condition = col_character(),
  date = col_date(format = "")
)
)

# Fit light curves for each individual
filmy_lc_model_fits <-
  filmy_lc %>%
  select(species, generation, individual, coll_num, condition, date, par, etr) %>%
  nest(data = c(par, etr)) %>%
  mutate(
    nls_mod = map(
      data, 
      # nonlinear least-squares estimates of parameters of the subsetted data
      ~nls(
        etr~max(etr)*(1-exp(-k*par)),
        start = list(k = 0.04),
        data =.,
        trace = FALSE,
        control = list(maxiter = 500))
    ),
    k_stats = map(nls_mod, broom::tidy)
  ) %>%
  mutate(fitted = map(nls_mod, ~fitted(.) %>% tibble(etr_fit = .))) %>%
  select(-nls_mod, -k_stats) %>%
  unnest(cols = c(data, fitted)) %>%
  # Convert genus to just first letter
  separate(species, c("genus", "epithet")) %>%
  mutate(genus = substr(genus, 1, 1)) %>%
  unite("species", genus, epithet) %>%
  mutate(id = paste3(species, coll_num, individual, as.character(date), sep = "_")) %>%
  # Nest fitted data
  group_by(species, generation, individual, coll_num, condition, date, id) %>%
  nest() %>%
  ungroup

# Subset to gametophytes
gameto_lc_model_fits <-
  filmy_lc_model_fits %>%
  filter(generation == "gametophyte")

# In RStudio, interactively go through data by increasing row selected one at a time,
# select outliers, and save each as CSV
select_lc_points(gameto_lc_model_fits, 5)
