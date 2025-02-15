# Clean light-curve raw data

# Load packages
library(readxl)
source("R/packages.R")

# Load functions
source("R/functions.R")

# Read in specimen collection data
# (requires clean_specimen_data.R to be run first)
specimens <- read_csv("data/fern_specimens.csv") 

# Read in "official" filmy species names (from growth habit data)
filmy_names <- read_csv("data/filmy_growth_habit.csv") %>% 
  select(species) %>%
  separate(species, c("genus", "specific_epithet"), remove = FALSE)

# Filmy fern sporophytes (lab measurements) ----

# Clean voucher data for filmy fern light curves measured in the lab
# column "sporo_site" provides name of site where the sporophyte sample
# was collected, if no voucher was made.
filmy_lc_lab_voucher_data <-
  read_excel("data_raw/2013/Filmy Fern Light Curves 7_26_13.xlsx", skip = 1) %>%
  clean_names %>%
  select(genus, specific_epithet = species, notes) %>%
  remove_empty("rows") %>%
  mutate(
    sporo_site = notes,
    sporo_site = str_remove_all(sporo_site, ", Nitta 3090"),
    notes = replace_na(notes, "none"),
    genus = ifelse(specific_epithet == "Hymenophyllum sp (3171)", "Hymenophyllum", genus),
    specific_epithet = str_replace_all(specific_epithet, "3178 A", "3178A"),
    specific_epithet = ifelse(specific_epithet == "Hymenophyllum sp (3171)", "sp (3171)", specific_epithet),
    # Crepidomanes minutum from Aorai 1700m Nitta '3090' was in error: should be 3190
    specific_epithet = ifelse(str_detect(notes, "Nitta 3090"), glue("{specific_epithet} (3190)"), specific_epithet),
    specific_epithet = str_remove_all(specific_epithet, "cf |-")
  ) %>%
  filter(!is.na(genus)) %>%
  separate(specific_epithet, c("specific_epithet", "coll_num"), extra = "drop", fill = "right") %>%
  select(genus, specific_epithet, coll_num, sporo_site) %>%
  mutate(species_manual = glue("{genus}_{specific_epithet}")) %>%
  left_join(transmute(specimens, species_from_coll_data = str_replace_all(species, " ", "_"), coll_num), by = "coll_num") %>%
  mutate(species = case_when(
    is.na(species_from_coll_data) ~ species_manual,
    TRUE ~ species_from_coll_data
  )) %>%
  select(species, coll_num, sporo_site)

filmy_sporo_lc_files <- c(
  "data_raw/2013/minipam/Filmy Fern  light curves 8-1-13.pam",
  "data_raw/2013/minipam/Filmy Fern  light curves 8-6-13.pam",
  "data_raw/2013/minipam/Filmy Fern  light curves 8-8-13.pam",
  "data_raw/2013/minipam/Abrodictyum_dentatum light curves 7-31-13.pam"
)

filmy_sporo_lc_lab <- 
  filmy_sporo_lc_files %>%
  set_names(fs::path_file(.)) %>%
  map(~parse_pam(., ret_type = "lc")) %>% 
  imap(~mutate(.x, file = .y)) %>%
  bind_rows() %>%
  # Filter out non-labeled data, standardize names
  mutate(num_only = str_detect(id, "[a-z]|[A-Z]", negate = TRUE)) %>%
  filter(!is.na(id), !num_only) %>%
  select(-num_only) %>%
  mutate(
    id = id %>%
      str_remove_all("\\.| l-int3") %>%
      str_replace_all("Abrodoctyum_cfasaegrayi", "Abrodictyum_asaegrayi1") %>%
      str_replace_all("C_apiifolia", "Callistopteris_apiifolia") %>%
      str_replace_all("H_multifidum", "Hymenophyllum_multifidum") %>%
      str_replace_all("H_palladum", "Hymenophyllum_pallidum") %>%
      str_replace_all("Hymenophyllum_sp_3171", "Hymenophyllum_braithwaitei") %>%
      str_replace_all("V_maxima", "Vandenboschia_maxima") %>%
      # Crepidomanes_minutum_3090 is in error for Crepidomanes_minutum_3190
      str_replace_all("Crepidomanes_minutum_3090", "Crepidomanes_minutum3") %>%
      str_replace_all("Crepidomanes_minutum_", "Crepidomanes_minutum2_")
    ) %>% 
  separate(id, c("genus", "epithet", "individual"), remove = FALSE, fill = "right") %>%
  unite("species", genus:epithet) %>%
  select(-id) %>%
  mutate(generation = "sporophyte") %>%
  # exclude an extra light curve for C. humile
  filter(!(species == "Crepidomanes_humile" & file == "Filmy Fern  light curves 8-1-13.pam")) %>%
  # exclude D. tahitense, which only had 1 replicate
  filter(species != "Didymoglossum_tahitense") %>%
  # add voucher data 
  left_join(filmy_lc_lab_voucher_data, by = "species") %>%
  select(-file) %>%
  mutate(condition = "lab") %>%
  # time stamp should always be unique
  assert(is_uniq, date_time)

# Field measurements ----

# Some of the rows in the .pam files have errors (e.g., light curve started
# but didn't finish), so manually inspect and exclude these
# so that the data can be read in properly.
field_lc_files <- tribble(
  ~file, ~exclude_lines,
  "data_raw/2013/minipam/Atiati 400m 7-4-13.pam", c(360:363, 216:225),
  "data_raw/2013/minipam/Mouaputa 200m 6-8-13.pam", NULL,
  "data_raw/2013/minipam/Mouaputa 400m 6-7-13.pam", c(123:126, 248:252, 253:262),
  "data_raw/2013/minipam/Mouaputa 400m 7-1-13.pam", NULL,
  "data_raw/2013/minipam/Mouaputa 600m 6-12-13.pam", NULL,
  "data_raw/2013/minipam/Mouaputa 600m 6-24-13.pam", NULL,
  "data_raw/2013/minipam/Mouaputa 800m 6-10-13.pam", c(190:200),
  "data_raw/2013/minipam/Mouaputa 800m 6-21-13.pam", NULL,
  "data_raw/2013/minipam/Mouaputa 800m 7-10-13.pam", NULL,
  "data_raw/2013/minipam/Muaroa 400m 6-19-13.pam", c(150:161, 346:349),
  "data_raw/2013/minipam/Muaroa 400m 6-23-13.pam", c(92:102, 104:116),
  "data_raw/2013/minipam/Rotui 600m 7-6-13.pam", c(58:69, 227:228),
  "data_raw/2013/minipam/Rotui 800m 7-8-13.pam", NULL,
  "data_raw/2013/minipam/Three Pines 6-30-13.pam", c(59:70, 214:227),
  "data_raw/2013/minipam/Three Pines 200m 6-16-13.pam", c(71:82, 328:333),
  "data_raw/2013/minipam/Three Pines 200m 6-28-13.pam", c(294:296, 384:385),
  "data_raw/2013/minipam/Three Pines 200m 6-29-13.pam", NULL
) %>%
  mutate(ret_type = "lc")

# Loop over the data files, extract light curves, then
# combine into a single dataframe including sporo_site and date
# extracted from file name.
# This includes some duplicates (in Three Pines 200m 6-28-13.pam and Three Pines 200m 6-29-13.pam)
field_lc_data_with_dups <-
field_lc_files %>%
  mutate(lc_data = pmap(., parse_pam)) %>%
  mutate(
    file_short = fs::path_file(file),
    sporo_site = str_replace_all(file_short, "Three Pines 6-30-13", "Three Pines 200m 6-30-13") %>%
      str_match("(^.*m) ") %>% 
      magrittr::extract(,2)) %>%
  select(lc_data, sporo_site) %>%
  unnest(cols = c(lc_data)) %>%
  mutate(condition = "field") %>%
  unique()

# Pull out the duplicates, and make them unique
de_duplicated_samples <-
  field_lc_data_with_dups %>%
  filter(!is_uniq(date_time)) %>%
  arrange(date_time, id) %>%
  mutate(is_dup = duplicated(date_time)) %>%
  filter(!is_dup) %>%
  select(-is_dup) %>%
  assert(is_uniq, date_time)
  
# Remove any duplicated data, and add back in the de-duplicated version
field_lc_data <-
field_lc_data_with_dups %>%
  filter(is_uniq(date_time)) %>%
  bind_rows(de_duplicated_samples) %>%
  # time stamp should always be unique
  assert(is_uniq, date_time)

# Filter to filmy fern sporophytes
filmy_genera <- filmy_names %>% pull(genus) %>% unique()

filmy_sporo_lc_field <-
field_lc_data %>%
  filter(str_detect(id, paste0(filmy_genera, collapse =  "|"))) %>%
  filter(!str_detect(id, regex("fail", ignore_case = TRUE))) %>%
  filter(!str_detect(id, regex("gameto", ignore_case = TRUE))) %>%
  # Exclude Crepidomanes minutum, because didn't record variety and have no way to verify
  filter(!str_detect(id, "minutum")) %>%
  mutate(id = str_remove_all(id, "2422_|_sporo|2507_")) %>%
  separate(id, c("genus", "epithet", "individual"), fill = "right", extra = "drop") %>%
  assert(in_set(filmy_genera), genus) %>%
  unite("species", genus, epithet) %>%
  mutate(generation = "sporophyte")
  
# Filter to filmy fern gametophytes
filmy_gameto_lc_field <-
  field_lc_data %>%
  # If the id code is less than 3 digits, it's not a collection number. Filter these out.
  mutate(coll_num = str_split(id, "_") %>% map_chr(1)) %>%
  filter(nchar(coll_num) > 2) %>% 
  # Filter out failures
  filter(!(str_detect(id, regex("fail", ignore_case = TRUE)))) %>% 
  filter(!is.na(yield)) %>%
  # Join with specimen data by collection number
  left_join(select(specimens, species, family, coll_num, generation), by = "coll_num") %>%
  filter(family == "Hymenophyllaceae") %>%
  filter(generation == "gametophyte") %>%
  rename(individual = id) %>%
  select(-family) %>%
  # Convert "individual" to a number code that is unique within each combination 
  # of species + collection number
  assert(not_na, coll_num, species) %>%
  group_by(coll_num) %>%
  mutate(individual = factor(individual) %>% numberify) %>%
  ungroup %>%
  mutate(individual = as.character(individual))

# Combine data ----
filmy_lc_data <- bind_rows(filmy_sporo_lc_lab, filmy_sporo_lc_field, filmy_gameto_lc_field) %>% 
  mutate(species = str_replace_all(species, " ", "_")) %>%
  # Check not missing data
  assert(not_na, type:generation, condition) %>%
  # Check species names formatted correctly
  assert(in_set(filmy_names$species), species) %>%
  # Check generation formatted correctly
  assert(in_set("sporophyte", "gametophyte"), generation) %>%
  # Check measurement condition formatted correctly
  assert(in_set("field", "lab"), condition) %>%
  # Check that no measurements are duplicated using time stamp
  assert(is_uniq, date_time) %>%
  # Add date
  mutate(date = date(date_time))

# Filter out poorly fitting light curves ----

# Fit model to each set of light curve data
filmy_lc_models <-
  filmy_lc_data %>%
  select(species, generation, individual, coll_num, condition, date, par, etr) %>%
  # create ID field
  mutate(id = paste(coll_num, individual, as.character(date), sep = "_")) %>%
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
  mutate(fitted = map(nls_mod, ~fitted(.) %>% tibble(etr_fit = .)))

# Extract model parameters: pvalue, critical PAR, and ETR at 95% of estimated max value
filmy_lc_model_params <-
  filmy_lc_models %>%
  select(-nls_mod, -data) %>%
  unnest(cols = fitted) %>%
  group_by(species, generation, individual, coll_num, condition, date) %>%
  mutate(etr_max = max(etr_fit)) %>%
  ungroup %>%
  select(-etr_fit) %>%
  unique() %>%
  unnest(cols = k_stats) %>%
  mutate(par_critical = -log(0.05)/estimate)

# Filter light curve data by removing any samples with a poor light curve fit
# (p-value > 0.05)
filmy_lc_data_filtered <-
  filmy_lc_data %>%
  # Add light-curve p-values
  left_join(
      select(filmy_lc_model_params, species, generation, individual, coll_num, condition, date, p.value)
  ) %>%
  # Make sure the join didn't change the number of rows
  verify(nrow(.) == nrow(filmy_lc_data)) %>%
  assert(not_na, p.value) %>%
  # Filter by p-value
  filter(p.value < 0.05)

filmy_lc_models_filtered <-
  filmy_lc_models %>%
  # Add light-curve p-values
  left_join(
    select(filmy_lc_model_params, species, generation, individual, coll_num, condition, date, p.value)
  ) %>%
  # Make sure the join didn't change the number of rows
  verify(nrow(.) == nrow(filmy_lc_models)) %>%
  # Filter by p-value
  assert(not_na, p.value) %>%
  filter(p.value < 0.05) %>%
  select(-p.value)

# Flag and filter outliers ----

# Optional: plot all filtered light curves to identify those with outliers
# - Extract fitted data points, abbreviate species names
filmy_lc_model_fitted_data_filtered <-
  filmy_lc_models_filtered %>%
  select(-nls_mod, -k_stats) %>%
  unnest(cols = c(data, fitted)) %>%
  abbrev_sp()

filmy_lc_model_params_filtered <-
  filmy_lc_model_params %>%
  filter(p.value < 0.05)  %>%
  abbrev_sp()

# - Optional: make plots to visualize raw light curves
filmy_all_lc_plot_filtered <-
  ggplot(abbrev_sp(filmy_lc_model_fitted_data_filtered), aes(x = par)) +
  geom_point(aes(y = etr)) +
  geom_line(aes(y = etr_fit)) +
  geom_hline(data = filmy_lc_model_params_filtered, aes(yintercept = etr_max), color = "green") +
  geom_vline(data = filmy_lc_model_params_filtered, aes(xintercept = par_critical), color = "red") +
  facet_wrap(vars(generation, species, id), scales = "free")

ggsave(plot = filmy_all_lc_plot_filtered, file = "results/all_light_curves_filtered.pdf", height = 40, width = 40)

# INTERACTIVE STEP:
# To run, uncomment the lines below
# In RStudio, go through data by increasing row selected one at a time,
# select outliers, and save each as CSV

# filmy_lc_models_filtered %>%
#   filter(str_detect(id, "2479_1_2013-06-10")) %>%
#   select_lc_points()
  
# Read in outliers, flag them in the filtered data

read_outlier <- function(file) {
  suppressWarnings(
    read_csv(file, col_types = cols(
      X1 = col_double(),
      species = col_character(),
      generation = col_character(),
      individual = col_character(),
      coll_num = col_character(),
      condition = col_character(),
      date = col_date(format = ""),
      par = col_double(),
      etr = col_double(),
      etr_fit = col_double()
    ))
  ) %>%
    select(-X1, -etr_fit)
}

outliers <- list.files("data_raw/intermediates/lc_outliers", full.names = TRUE) %>%
  map_df(read_outlier) %>%
  unique() %>%
  mutate(outlier = TRUE)

filmy_lc_data_filtered_outliers_flagged <-
  filmy_lc_data_filtered %>%
  left_join(
    outliers,
    by = c("par", "etr", "species", "individual", "generation", "coll_num", "condition", "date")) %>%
  # make sure join didn't duplicate any values
  # (number of rows should be the same before and after)
  verify(nrow(.) == nrow(filmy_lc_data_filtered)) %>%
  mutate(outlier = replace_na(outlier, FALSE)) %>%
  # Add a light-curve ID number for each sample:
  # unique to combination of: species, generation, individual, coll_num, condition, date
  # - replace NA values in collection number with "sn" (no number)
  mutate(coll_num = replace_na(coll_num, "sn")) %>%
  assert(not_na, species, generation, individual, coll_num, condition, date, yield, par, etr, outlier) %>%
  group_by(species, generation, individual, coll_num, condition, date) %>%
  mutate(light_id = cur_group_id()) %>%
  ungroup() %>%
  assert(not_na, light_id) %>%
  # Don't need p-value any more after flagging outliers
  select(-p.value)

# Write out final cleaned dataset ----
write_csv(filmy_lc_data_filtered_outliers_flagged, "data/filmy_light_curves.csv")
