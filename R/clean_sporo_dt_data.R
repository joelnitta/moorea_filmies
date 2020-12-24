# Clean sporophyte desiccation tolerance (DT) test raw data

# Load packages
library(readxl)
source("R/packages.R")

# Load functions
source("R/functions.R")

# 2012 data ----
### Fix column names ###
# In the original excel file, these are in two rows.
moorea_filmy_dt_2012_raw_path <- "data_raw/2012/Moorea Filmy DT 2 8-26 (version 1).xlsx"

headers_1 <- read_excel(moorea_filmy_dt_2012_raw_path) %>% colnames()

headers_2 <- read_excel(moorea_filmy_dt_2012_raw_path, skip = 1) %>% colnames()

dt_2012_head_1 <- read_excel(moorea_filmy_dt_2012_raw_path, skip = 2, col_names = headers_1) %>%
  remove_empty(which = "cols", quiet = FALSE) %>% # delete empty columns
  clean_names()

dt_2012_head_2 <- read_excel(moorea_filmy_dt_2012_raw_path, skip = 2, col_names = headers_2) %>%
  remove_empty(which = "cols", quiet = FALSE) %>%
  clean_names()

headers <- tibble(
  header_1 = colnames(dt_2012_head_1),
  header_2 = colnames(dt_2012_head_2)
)

# Write out these column headers, then manually edit new column names
# write.csv(headers, "moorea_filmy_dt_2012_headers.csv")

# Read in file with new column names
new_headers <- read_csv("data_raw/intermediates/moorea_filmy_dt_2012_headers.csv")

dt_2012_all_raw <- read_excel(moorea_filmy_dt_2012_raw_path, skip = 2, col_names = FALSE) %>%
  remove_empty(which = "cols", quiet = FALSE) %>%
  set_names(new_headers$new_header)

### Fix contents ###
dt_2012_all <-
  dt_2012_all_raw %>% 
  select(
    species,
    salt,
    dry_time,
    individual,
    matches("weight|yield")
  ) %>%
  # Fix species names
  mutate(
    species = case_when(
    species == "C bip" ~ "Crepidomanes_bipunctatum",
    species == "C hum" ~ "Crepidomanes_humile",
    species == "C min" ~ "Crepidomanes_minutum1",
    species == "Cal api" ~ "Callistopteris_apiifolia",
    species == "H pol" ~ "Hymenophyllum_polyanthos",
    species == "P end" ~ "Polyphlebium_borbonicum",
    TRUE ~ NA_character_
  )) %>%
  # Convert dry time to factor (control will be NA)
  mutate(
    dry_time = case_when(
      dry_time == "15" ~ 15,
      dry_time == "2" ~ 2,
      TRUE ~ NaN
    ),
    dry_time = as.factor(dry_time)
  ) %>%
  mutate(salt = str_replace_all(salt, "control", "Control")) %>%
  # Weight and light measures should be numeric. 
  # Some of these had string comments (e.g., "missing") if the value was missing.
  # parse_number() will turn those comments to NA.
  mutate(across(matches("yield|weight") & where(is.character), parse_number)) %>%
  mutate(weight_desiccated_2d = ifelse(dry_time == 2, weight_desiccated, NA)) %>%
  mutate(weight_desiccated_15d = ifelse(dry_time == 15, weight_desiccated, NA)) %>%
  select(-weight_desiccated)

# Split into test and control data to reshape further
dt_2012_test <-
dt_2012_all %>%
  filter(salt != "Control") %>%
  remove_empty(which = "cols", quiet = TRUE)

# The same individuals were used as control for both 2d and 15d.
# But when combining data, enter each into a separate row.
dt_2012_2d_control <-
dt_2012_all %>%
  filter(salt == "Control") %>%
  remove_empty(which = "cols", quiet = TRUE) %>% 
  select(
    species, salt, dry_time, individual, weight_pre, 
    weight_desiccated_2d = weight_desiccated_control_2d, 
    yield_pre, yield_30min, yield_24hr, yield_48hr) %>%
  mutate(dry_time = 2 %>% factor(levels = c(2, 15)))

dt_2012_15d_control <-
  dt_2012_all %>%
  filter(salt == "Control") %>%
  remove_empty(which = "cols", quiet = TRUE) %>% 
  select(
    species, salt, dry_time, individual, weight_pre, 
    weight_desiccated_15d = weight_desiccated_control_15d, 
    yield_pre, 
    yield_30min = yield_30min_control_15d, 
    yield_24hr = yield_24hr_control_15d,
    yield_48hr = yield_48hr_control_15d) %>%
  mutate(dry_time = 15 %>% factor(levels = c(2, 15)))

filmy_dt_2012 <- bind_rows(
  dt_2012_test,
  dt_2012_2d_control,
  dt_2012_15d_control
) %>%
  # combine weight_desiccated_2d, weight_desiccated_15d into a single column
  mutate(
    weight_desiccated = case_when(
      dry_time == "2" ~ weight_desiccated_2d,
      dry_time == "15" ~ weight_desiccated_15d
    )
  ) %>% 
  select(-c(weight_desiccated_2d, weight_desiccated_15d)) %>%
  mutate(dataset = "2012") %>%
  check_dt_data

# 2013 data ----

filmy_dt_2013_1_raw <- read_excel("data_raw/2013/Filmy Fern DT 1 Dtah, Ckur, Capi.xlsx", sheet = 1)

filmy_dt_2013_1 <- filmy_dt_2013_1_raw %>%
  clean_names() %>%
  # remove a single Callistopteris gametophyte control sample
  filter(!str_detect(species, "gametophyte")) %>%
  remove_empty(which = "rows", quiet = TRUE) %>%
  transmute(
    species,
    salt = salt_treatment, 
    dry_time = parse_number(duration_of_treatment),
    individual,
    yield_pre = wet_yield,
    memory_pre = wet_last_memory,
    yield_30min = x30min_yield,
    memory_30min = x30min_last_memory,
    yield_24hr = x24hr_mini_pam_data,
    memory_24hr = x24_hr_last_memory,
    yield_48hr = x48hr_mini_pam_data,
    memory_48hr = x48_hr_last_memory,
    yield_72hr = x72hr_mini_pam_c_apii_only,
    memory_72hr = x72_hr_last_memory) %>% 
  # Fix species names
  mutate(
    species = case_when(
      species == "D" ~ "Didymoglossum_tahitense",
      species == "Cr" ~ "Crepidomanes_kurzii",
      species == "Cl" ~ "Callistopteris_apiifolia",
      species == "Didymoglossum" ~ "Didymoglossum_tahitense",
      species == "Crepidomanes kurzii" ~ "Crepidomanes_kurzii",
      TRUE ~ NA_character_
    )) %>%
  # Convert dry time to factor (control will be NA)
  mutate(
    dry_time = case_when(
      dry_time == "15" ~ 15,
      dry_time == "2" ~ 2,
      TRUE ~ NaN
    ),
    dry_time = as.factor(dry_time)
  ) %>%
  mutate(
    yield_pre = ifelse(salt != "Control", 1000*yield_pre, yield_pre),
    dataset = "2013_1",
    across(contains("memory") & where(is.character), parse_number)) %>%
  check_dt_data

filmy_dt_2013_2_raw <- read_excel("data_raw/2013/Filmy Fern DT 2 Pbor.xlsx", sheet = 1)

filmy_dt_2013_2 <- filmy_dt_2013_2_raw %>%
  clean_names() %>%
  remove_empty(which = "rows", quiet = TRUE) %>%
  transmute(
    species,
    salt = salt_treatment, 
    dry_time = duration_of_treatment,
    individual,
    yield_pre = yield_wet,
    memory_pre = memory_wet,
    yield_30min = yield_20_min_recover,
    memory_30min = memory_20_min_recover,
    yield_24hr = yield_24hr_recover,
    memory_24hr = memory_24hr_recover,
    yield_48hr = yield_48hr_recover,
    memory_48hr = memory_48hr_recover,
    weight_pre = wet_weight_g,
    weight_desiccated_2d = desiccated_weight_2d_g,
    weight_desiccated_15d = desiccated_weight_15d,
    weight_dry = dry_weight_g) %>% 
  # Fix species names
  mutate(
    species = case_when(
      species == "P. borbonicum" ~ "Polyphlebium_borbonicum",
      TRUE ~ NA_character_
    )) %>%
  assert(not_na, species) %>% 
  # Convert dry time to factor
  mutate(dry_time = as.factor(dry_time)) %>%
  # combine weight_desiccated_2d, weight_desiccated_15d into a single column
  mutate(
    weight_desiccated = case_when(
      dry_time == "2" ~ weight_desiccated_2d,
      dry_time == "15" ~ weight_desiccated_15d
    )
  ) %>% 
  select(-c(weight_desiccated_2d, weight_desiccated_15d)) %>%
  mutate(
    dataset = "2013_2",
    across(contains("memory") & where(is.character), parse_number)) %>%
  check_dt_data

filmy_dt_2013_3_raw <- read_excel("data_raw/2013/Filmy Fern DT 3 Hmul, Hdig.xlsx", sheet = 1)

filmy_dt_2013_3 <- filmy_dt_2013_3_raw %>%
  clean_names() %>%
  remove_empty(which = "rows", quiet = TRUE) %>%
  transmute(
    species,
    salt = salt_treatment, 
    dry_time = parse_number(duration_of_treatment),
    individual,
    yield_pre = mini_pam_data_wet,
    memory_pre = last_memory_for_wet_data,
    yield_30min = mini_pam_data_30min,
    memory_30min = last_memory_for_30_min_data,
    yield_24hr = mini_pam_data_24hr,
    memory_24hr = last_memory_for_24hr,
    yield_48hr = parse_number(mini_pam_data_48hr),
    memory_48hr = last_memory_for_48hr,
    weight_pre = wet_weight_g,
    weight_desiccated_2d = x2d_desiccated_weight_g,
    weight_desiccated_15d = dessicated_weight,
    weight_dry = oven_dry_weight_g) %>%
  # Fix species names
  mutate(
    species = case_when(
      species == "H. digitatum" ~ "Hymenophyllum_digitatum",
      species == "H. multifidum" ~ "Hymenophyllum_multifidum",
      TRUE ~ NA_character_
    )) %>%
  assert(not_na, species) %>%
  # Convert dry time to factor (control will be NA)
  mutate(dry_time = as.factor(dry_time)) %>%
  # Combine weight_desiccated_2d, weight_desiccated_15d into a single column
  mutate(
    weight_desiccated = case_when(
      dry_time == "2" ~ weight_desiccated_2d,
      dry_time == "15" ~ weight_desiccated_15d
    )
  ) %>% 
  select(-c(weight_desiccated_2d, weight_desiccated_15d)) %>%
  mutate(
    dataset = "2013_3",
    across(contains("memory") & where(is.character), parse_number)) %>%
  check_dt_data

filmy_dt_2013_4_raw <- read_excel("data_raw/2013/Filmy Fern DT 4 Hpol.xlsx", sheet = 1)

filmy_dt_2013_4 <- filmy_dt_2013_4_raw %>% 
  clean_names() %>%
  filter(desiccated_weight_g != "(48 hour weight)") %>%
  transmute(
    species,
    salt = salt_treatment, 
    dry_time = parse_number(duration_of_treatment),
    individual,
    yield_pre = mini_pam_data_wet,
    memory_pre = last_memory_for_wet_data,
    yield_30min = yield_20_min_recover,
    memory_30min = memory_20_min,
    yield_24hr = parse_number(yield_24hr_recover),
    memory_24hr = memory_24hr_recover,
    yield_48hr = yield_48hr_recovery,
    memory_48hr = memory_48hr_recovery,
    weight_pre = wet_weight_g,
    weight_desiccated_2d = parse_number(desiccated_weight_g),
    weight_desiccated_15d = desiccated_weight_15d,
    weight_dry = oven_dry_weight_g) %>%
  remove_empty(which = "rows", quiet = TRUE) %>%
  remove_empty(which = "cols", quiet = TRUE) %>%
  # Fix species names
  mutate(
    species = case_when(
      species == "H. Polyanthos" ~ "Hymenophyllum_polyanthos",
      TRUE ~ NA_character_
    )) %>% 
  # Convert dry time to factor (control will be NA)
  mutate(dry_time = as.factor(dry_time)) %>%
  # Combine weight_desiccated_2d, weight_desiccated_15d into a single column
  mutate(
    weight_desiccated = case_when(
      dry_time == "2" ~ weight_desiccated_2d,
      dry_time == "15" ~ weight_desiccated_15d
    )
  ) %>% 
  select(-c(weight_desiccated_2d, weight_desiccated_15d)) %>%
  mutate(
    dataset = "2013_4",
    across(contains("memory") & where(is.character), parse_number)) %>%
  check_dt_data

filmy_dt_2013_5_raw <- read_excel("data_raw/2013/Filmy Fern DT Set 5 DONE.xlsx", sheet = 1)

filmy_dt_2013_5 <- filmy_dt_2013_5_raw %>%
  clean_names() %>%
  transmute(
    species,
    salt = salt_treatment, 
    # some of these are blank, but A. dentatum was only tested for two days
    dry_time = replace_na(duration_of_treatment, 2),
    individual,
    yield_pre = mini_pam_data_wet,
    memory_pre = last_memory_for_wet_data_7_26,
    # fvfm_dry = yield_dry,
    yield_30min = yield_20_min_recover_7_28_18_50pm,
    memory_30min = memory_20_min,
    yield_24hr = yield_24hr_recover_7_29,
    memory_24hr = memory_24hr_recover,
    yield_48hr = yield_48hr_recovery_7_30,
    memory_48hr = memory_48hr_recovery,
    weight_pre = wet_weight_g,
    weight_desiccated = desiccated_weight_2d_g,
    weight_30min = mass_20min_recover,
    weight_24hr = weight_24hr,
    weight_48hr = weight_48_hr,
    weight_dry = oven_dry_weight_g) %>%
  remove_empty(which = "rows", quiet = TRUE) %>%
  remove_empty(which = "cols", quiet = TRUE) %>%
  # Fix species names
  mutate(
    species = case_when(
      species == "A dentatum" ~ "Abrodictyum_dentatum",
      TRUE ~ NA_character_
    )) %>% 
  # Convert dry time to factor (control will be NA)
  mutate(
    dry_time = as.factor(dry_time),
    dataset = "2013_5",
    across(contains("memory") & where(is.character), parse_number)
  ) %>%
  check_dt_data

# 2014 data ----

filmy_dt_2014_raw <- read_excel("data_raw/2014/filmyDT2014.xlsx", sheet = 1)

filmy_dt_2014 <- filmy_dt_2014_raw %>%
  clean_names() %>%
  transmute(
    species,
    salt, 
    # some of these are blank, but A. dentatum was only tested for two days
    dry_time = time,
    individual,
    yield_pre = wet_fv_fm,
    memory_pre = wet_mem,
    # fvfm_dry = dry_fv_fm,
    yield_30min = x30min_fv_fm,
    memory_30min = x30min_mem,
    yield_24hr = x24hr_fv_fm,
    memory_24hr = x24hr_mem,
    yield_48hr = x48hr_fv_fm,
    memory_48hr = x48hr_mem,
    weight_pre = wet_weight,
    weight_desiccated = dry_weight_7_2_14,
    weight_30min = x30min_weight_7_2_14,
    weight_24hr = x24hr_weight_7_2_14,
    weight_48hr = x48hr_weight_7_2_14,
    weight_dry = oven_dry_weight) %>%
  remove_empty(which = "rows", quiet = TRUE) %>%
  remove_empty(which = "cols", quiet = TRUE) %>%
  # Fix species names
  mutate(
    species = case_when(
      species == "C. minutum v2" ~ "Crepidomanes_minutum2",
      species == "H. pallidum" ~ "Hymenophyllum_pallidum",
      species == "V. maximum" ~ "Vandenboschia_maxima",
      TRUE ~ NA_character_
    )) %>% 
  # Convert dry time to factor (control will be NA)
  mutate(dry_time = dry_time %>% parse_number %>% as.factor,
         dataset = "2014",
         across(contains("memory") & where(is.character), parse_number)) %>%
  check_dt_data

# Combine data from different years ----
filmy_sporo_dt <- bind_rows(
  filmy_dt_2012, 
  filmy_dt_2013_1, 
  filmy_dt_2013_2, 
  filmy_dt_2013_3, 
  filmy_dt_2013_4, 
  filmy_dt_2013_5, 
  filmy_dt_2014) %>%
  select(
    species, salt, dry_time, individual, dataset, contains("weight"), contains("yield"), contains("memory")
  ) %>%
  mutate(generation = "sporophyte") %>%
  check_dt_data

# Join with PAM data to add timestamps ----
# Convert to long format
filmy_sporo_dt %>%
  select(-contains("yield")) %>%
  pivot_longer(names_to = "condition", values_to = "memory", contains("memory")) %>%
  filter(!is.na(memory)) %>%
  mutate(condition = str_remove_all(condition, "memory_") %>% 
           factor(levels = c("pre", "30min", "24hr", "48hr", "72hr"))
  ) %>%
  arrange(dry_time, condition, salt, species, individual) %>%
  mutate(pam_set = case_when(
    dry_time == 15 & condition %in% c("30min", "24hr", "48hr", "72hr") ~ "set2",
    TRUE ~ "set1"
  )) %>%
  assert_rows(col_concat, is_uniq, pam_set, memory)


# Write to data_raw/ ----
write_csv(filmy_sporo_dt, "data/filmy_sporo_dt.csv")
