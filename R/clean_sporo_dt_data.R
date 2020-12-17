# Load packages
source(here::here("R/packages.R"))

# Custom funcs ---
source(here("R/functions.R"))

# 2012 data ----
### Fix column names ###
# In the original excel file, these are in two rows.
moorea_filmy_dt_2012_raw_path <- "data/2012/Moorea Filmy DT 2 8-26 (version 1).xlsx"

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

# write.csv(headers, "moorea_filmy_dt_2012_headers.csv")

new_headers <- read_csv("data/intermediates/moorea_filmy_dt_2012_headers.csv")

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

filmy_dt_2013_1_raw <- read_excel("data/2013/Filmy Fern DT 1 Dtah, Ckur, Capi.xlsx", sheet = 1)

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
    yield_30min = x30min_yield,
    yield_24hr = x24hr_mini_pam_data,
    yield_48hr = x48hr_mini_pam_data,
    yield_72hr = x72hr_mini_pam_c_apii_only) %>% 
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
  mutate(yield_pre = ifelse(salt != "Control", 1000*yield_pre, yield_pre)) %>%
  mutate(dataset = "2013_1") %>%
  check_dt_data

filmy_dt_2013_2_raw <- read_excel("data/2013/Filmy Fern DT 2 Pbor.xlsx", sheet = 1)

filmy_dt_2013_2 <- filmy_dt_2013_2_raw %>%
  clean_names() %>%
  remove_empty(which = "rows", quiet = TRUE) %>%
  transmute(
    species,
    salt = salt_treatment, 
    dry_time = duration_of_treatment,
    individual,
    yield_pre = yield_wet,
    yield_30min = yield_20_min_recover,
    yield_24hr = yield_24hr_recover,
    yield_48hr = yield_48hr_recover,
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
  mutate(dataset = "2013_2") %>%
  check_dt_data

filmy_dt_2013_3_raw <- read_excel("data/2013/Filmy Fern DT 3 Hmul, Hdig.xlsx", sheet = 1)

filmy_dt_2013_3 <- filmy_dt_2013_3_raw %>%
  clean_names() %>%
  remove_empty(which = "rows", quiet = TRUE) %>%
  transmute(
    species,
    salt = salt_treatment, 
    dry_time = parse_number(duration_of_treatment),
    individual,
    yield_pre = mini_pam_data_wet,
    yield_30min = mini_pam_data_30min,
    yield_24hr = mini_pam_data_24hr,
    yield_48hr = parse_number(mini_pam_data_48hr),
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
  mutate(dataset = "2013_3") %>%
  check_dt_data

filmy_dt_2013_4_raw <- read_excel("data/2013/Filmy Fern DT 4 Hpol.xlsx", sheet = 1)

filmy_dt_2013_4 <- filmy_dt_2013_4_raw %>% 
  clean_names() %>%
  filter(desiccated_weight_g != "(48 hour weight)") %>%
  transmute(
    species,
    salt = salt_treatment, 
    dry_time = parse_number(duration_of_treatment),
    individual,
    yield_pre = mini_pam_data_wet,
    yield_30min = yield_20_min_recover,
    yield_24hr = parse_number(yield_24hr_recover),
    yield_48hr = yield_48hr_recovery,
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
  mutate(dataset = "2013_4") %>%
  check_dt_data

filmy_dt_2013_5_raw <- read_excel("data/2013/Filmy Fern DT Set 5 DONE.xlsx", sheet = 1)

filmy_dt_2013_5 <- filmy_dt_2013_5_raw %>%
  clean_names() %>% 
  transmute(
    species,
    salt = salt_treatment, 
    # some of these are blank, but A. dentatum was only tested for two days
    dry_time = replace_na(duration_of_treatment, 2),
    individual,
    yield_pre = mini_pam_data_wet,
    # fvfm_dry = yield_dry,
    yield_30min = yield_20_min_recover_7_28_18_50pm,
    yield_24hr = yield_24hr_recover_7_29,
    yield_48hr = yield_48hr_recovery_7_30,
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
  mutate(dry_time = as.factor(dry_time)) %>%
  mutate(dataset = "2013_5") %>%
  check_dt_data

# 2014 data ----

filmy_dt_2014_raw <- read_excel("data/2014/filmyDT2014.xlsx", sheet = 1)

filmy_dt_2014 <- filmy_dt_2014_raw %>%
  clean_names() %>%
  transmute(
    species,
    salt, 
    # some of these are blank, but A. dentatum was only tested for two days
    dry_time = time,
    individual,
    yield_pre = wet_fv_fm,
    # fvfm_dry = dry_fv_fm,
    yield_30min = x30min_fv_fm,
    yield_24hr = x24hr_fv_fm,
    yield_48hr = x48hr_fv_fm,
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
  mutate(dry_time = dry_time %>% parse_number %>% as.factor) %>%
  mutate(dataset = "2014") %>%
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
    species, salt, dry_time, individual, dataset, contains("weight"), contains("yield")
  ) %>%
  mutate(generation == "sporo") %>%
  check_dt_data

### DT recovery ###
filmy_dt_recovery <-
  filmy_dt_wide %>%
  select(species, salt, dry_time, individual, matches("yield")) %>%
  select(-yield_72hr) %>%
  pivot_longer(names_to = "rec_time", values_to = "yield_recover", matches("30|24|48")) %>%
  mutate(
    rec_time = str_remove_all(rec_time, "yield_"),
    recovery = yield_recover / yield_pre)

### Weights ###
no_weight_data <- 
  filmy_dt_wide %>%
  select(species, salt, dry_time, individual, matches("weight")) %>%
  filter(across(matches("weight"), is.na))

filmy_dt_weight <-
  filmy_dt_wide %>%
  select(species, salt, dry_time, individual, matches("weight")) %>%
  anti_join(no_weight_data, by = c("species", "salt", "dry_time", "individual")) %>%
  mutate(percent_loss_2d =  1 - ((desiccated_weight_2d - oven_dry_weight) / (wet_weight - oven_dry_weight))) %>%
  mutate(percent_loss_15d =  1 - ((desiccated_weight_15d - oven_dry_weight) / (wet_weight - oven_dry_weight)))

# Visualize the data ----
filmy_dt_recovery_mean <-
filmy_dt_recovery %>%
  mutate(rec_time = factor(rec_time, levels = c("30min", "24hr", "48hr"), ordered = TRUE)) %>%
  group_by(species, salt, dry_time, rec_time) %>%
  summarize(
    recovery = mean(recovery, na.rm = FALSE),
    n = n(),
    .groups = "drop")

ggplot(filmy_dt_recovery_mean, aes(x = rec_time, y = recovery, group = interaction(dry_time, salt))) +
  geom_point(aes(shape = salt, color = salt)) +
  geom_line(aes(linetype = dry_time)) +
  facet_wrap(vars(species))