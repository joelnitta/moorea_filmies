# Clean gametophyte desiccation tolerance (DT) test raw data

# Load packages
library(readxl)
source("R/packages.R")

# Load functions
source("R/functions.R")

# Read in mini-PAM data ----

# Loop over the files: read them in, process each data frame in the list, then combine
minipam_2013_long <-
  # List all the 2013 gameto mini pam files
  list.files("data_raw/2013/minipam", pattern = "miniPAM SET\\d ALL.*\\.pam", full.names = TRUE) %>%
  set_names(fs::path_file(.)) %>%
  # Parse each one in a list
  map(
    ~parse_pam(., ret_type = "fl", recalc_yield = TRUE) %>%
      select(memory, yield, date_time, yield_error)) %>%
  # Add a column for "file"
  imap(~mutate(.x, file = .y)) %>%
  # Combine into single df
  bind_rows %>%
  # Add column for "pam set" (sets 1-5)
  mutate(pam_set = str_match(file, "SET(\\d)") %>% magrittr::extract(,2) %>% paste0("set", .)) %>%
  assert(not_na, pam_set) %>%
  select(-file) %>%
  filter(
    # Exclude rows with corrupted times (all should be 2013)
    date_time > lubridate::ymd_hms("2013-01-01 01:00:00"),
    date_time < lubridate::ymd_hms("2014-01-01 01:00:00")) %>%
  # Exclude some corrupted data in set 1
  filter(!(memory %in% c(1963, 1969:2001) & pam_set == "set1")) %>%
  arrange(date_time)

# 2012 sporophyte DT data ----
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

# 2013 sporophyte DT data ----

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
    across(contains("memory") & where(is.character), parse_number),
    # In a subset of yields, value was recorded as a decimal instead of whole number
    yield_30min = case_when(
      dry_time == "2" & salt != "H2O" ~ yield_30min*1000,
      TRUE ~ yield_30min
    ) ) %>%
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
    across(contains("memory") & where(is.character), parse_number))

# P. borbonicum 15d control was missing numbers for individual. Fill these in.
filmy_dt_2013_2$individual[filmy_dt_2013_2$"individual" %>% is.na() %>% which()] <- 1:8

check_dt_data(filmy_dt_2013_2)

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

# 2014 sporophyte DT data ----

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

# Combine sporophyte data from different years ----
filmy_sporo_dt_raw_yield <- bind_rows(
  filmy_dt_2012, 
  filmy_dt_2013_1, 
  filmy_dt_2013_2, 
  filmy_dt_2013_3, 
  filmy_dt_2013_4, 
  filmy_dt_2013_5, 
  filmy_dt_2014) %>%
  select(
    species, salt, dry_time, individual, dataset, matches("yield|memory|weight")
  ) %>%
  mutate(
    generation = "sporophyte",
    individual = as.character(individual)) %>%
  check_dt_data

# 2012 gametophyte DT data ----

# Read in each excel sheet as a separate dataframe.
# In 2012 only, there were three replicate measurements taken of each individual
# per treatment
gameto_2012_dark <- read_excel("data_raw/2012/all gametos 10-3.xlsx", sheet = "DARK", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x9, yield_3 = x13)

gameto_2012_dry <- read_excel("data_raw/2012/all gametos 10-3.xlsx", sheet = "DRY", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

gameto_2012_30min <- read_excel("data_raw/2012/all gametos 10-3.xlsx", sheet = "10 MIN", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

gameto_2012_24hr <- read_excel("data_raw/2012/all gametos 10-3.xlsx", sheet = "24 HR", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

gameto_2012_48hr <- read_excel("data_raw/2012/all gametos 10-3.xlsx", sheet = "48 HR", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

gameto_2012_72hr <- read_excel("data_raw/2012/all gametos 10-3.xlsx", sheet = "72 HR", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

# Do the same for an additional set of data from the 200m terrestrial site

gameto_2012_200m_ter_dark <- read_excel("data_raw/2012/new 200 m ter gameto.xlsx", sheet = "Pre", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x7, yield_2 = x12, yield_3 = x17)

gameto_2012_200m_ter_dry <- read_excel("data_raw/2012/new 200 m ter gameto.xlsx", sheet = "Dry + 20 min", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x11, yield_2 = x16, yield_3 = x21)

gameto_2012_200m_ter_30min <- read_excel("data_raw/2012/new 200 m ter gameto.xlsx", sheet = "Dry + 20 min", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x11, yield_2 = x10, yield_3 = x15)

gameto_2012_200m_ter_24hr <- read_excel("data_raw/2012/new 200 m ter gameto.xlsx", sheet = "24 hr", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x7, yield_2 = x12, yield_3 = x17)

gameto_2012_200m_ter_48hr <- read_excel("data_raw/2012/new 200 m ter gameto.xlsx", sheet = "48 hr", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x7, yield_2 = x12, yield_3 = x17)

gameto_2012_200m_ter_72hr <- read_excel("data_raw/2012/new 200 m ter gameto.xlsx", sheet = "72 hr", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x4, yield_2 = x9, yield_3 = x14)

# Combine the separate datasheets into a single wide dataframe including
# all replicate measurements
gameto_dt_2012_raw <- reduce(
  list(
  gameto_2012_dark %>% rename_with(~paste0(., "_pre"), yield_1:yield_3),
  gameto_2012_dry %>% rename_with(~paste0(., "_dry"), yield_1:yield_3),
  gameto_2012_24hr %>% rename_with(~paste0(., "_24hr"), yield_1:yield_3),
  gameto_2012_48hr %>% rename_with(~paste0(., "_48hr"), yield_1:yield_3),
  gameto_2012_72hr %>% rename_with(~paste0(., "_72hr"), yield_1:yield_3)
  ),
  left_join,
  by = c("plot", "individual")
) %>%
  select(-plot)

gameto_2012_200m_ter_raw <- reduce(
  list(
    gameto_2012_200m_ter_dark %>% rename_with(~paste0(., "_pre"), yield_1:yield_3),
    gameto_2012_200m_ter_dry %>% rename_with(~paste0(., "_dry"), yield_1:yield_3),
    gameto_2012_200m_ter_24hr %>% rename_with(~paste0(., "_24hr"), yield_1:yield_3),
    gameto_2012_200m_ter_48hr %>% rename_with(~paste0(., "_48hr"), yield_1:yield_3),
    gameto_2012_200m_ter_72hr %>% rename_with(~paste0(., "_72hr"), yield_1:yield_3)
  ),
  left_join,
  by = "individual"
)

# Obtain mean values for each individual/treatment
gameto_dt_2012 <-
  bind_rows(gameto_dt_2012_raw, gameto_2012_200m_ter_raw) %>%
  assert(is_uniq, individual) %>%
  pivot_longer(names_to = "var", values_to = "value", -individual) %>%
  separate(var, c("meas_type", "replicate", "condition"), sep = "_") %>%
  select(-meas_type) %>%
  group_by(individual, condition) %>%
  summarize(yield = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = "condition", values_from = "yield") %>%
  rename_with(~paste0("yield_", .), -individual) %>%
  mutate(salt = "MgNO3", dry_time = 2, individual = as.character(individual))

# 2013 gametophyte DT data ----

### Read in gametophyte DT test data ###

# These were originally entered on google drive, and have been
# downloaded from there as xlsx files.
# They have slightly differing column names, but the important
# ones can be identified by containing "yield" or "mem" (memory)

# Read in all gametophyte DT data as a list,
# but don't fix column names yet
gameto_dt_2013_raw_names <- 
  # List all 2013 gametophyte DT data sets 2-13 files
  list.files("data_raw/2013", pattern = regex("moorea 2013 gameto"), ignore.case = TRUE, full.names = TRUE) %>%
  # Name by vector of files by the file name, so we can include it when combining the data
  set_names(fs::path_file(.)) %>%
  # Read in the xlsx files as a list
  map(suppressMessages(read_excel)) %>%
  map(clean_names)

# Manually fix some column names that fix_gameto_dt_names() wouldn't handle properly
gameto_dt_2013_raw_names[["Moorea 2013 gametophyte DT set1.xlsx"]] <-
  gameto_dt_2013_raw_names[["Moorea 2013 gametophyte DT set1.xlsx"]] %>%
  select(
    individual = gametophytes,
    memory_pre = wet_memory,
    memory_dry = dry_last_memory,
    memory_30min = x20_min_mem,
    memory_24hr = x24_hr_mem,
    memory_48hr = x48_hr_mem,
    memory_72hr = x72_hr_mem
  ) %>%
  remove_empty("rows") %>%
  mutate(individual = make.unique(as.character(individual)))

gameto_dt_2013_raw_names[["Moorea 2013 Gameto DT Set 13 DONE.xlsx"]] <-
  gameto_dt_2013_raw_names[["Moorea 2013 Gameto DT Set 13 DONE.xlsx"]] %>%
  select(
    individual = gametophytes,
    memory_pre = memory,
    yield_pre = yield_wet_7_20,
    memory_dry = last_memory_in_dry_state,
    yield_dry = yield_in_dry_state_7_28,
    memory_30min = last_memory_rewetted,
    yield_30min = yield_rewetted_7_28,
    memory_24hr = x24_hr_mem,
    yield_24hr = x24_hr_yield_7_29,
    memory_48hr = x48_hr_mem,
    yield_48hr = x48_hr_yield_7_30,
    memory_72hr = x72_hr_mem,
    yield_72hr = x72_hr_yield_7_31
  )

# Fix column names and combine data
# Note that some data are missing in the raw data files, 
# not because fix_gameto_dt_names() didn't work properly:
# DT set 7, 8 is missing memory for pre-treatment
# DT set 6 is missing memory for dry state
# DT set 3 is missing memory and yield for dry state
gameto_dt_2013_raw_yield <- 
  gameto_dt_2013_raw_names %>%
  map(fix_gameto_dt_names) %>%
  # Add a column for "file"
  imap(~mutate(.x, file = .y)) %>%
  bind_rows() %>%
  mutate(
    # For some reason, a ".0" got appended to some collection numbers, remove this
    individual = str_remove_all(individual, "\\.0"),
    notes = replace_na(notes, "none")) %>%
  # Remove a row that is almost all empty except one erroneously entered value
  filter(!(is.na(individual) & str_detect(file, "Set 5"))) %>%
  # Remove rows with duplicated individuals, but differing measurements: not
  # sure which measurement is correct
  filter(!(individual %in% c("3105", "3106") & str_detect(file, "Set 12"))) %>%
  # Remove samples that were "missing"
  filter(
    !str_detect(notes, regex("missing", ignore_case = TRUE)),
    # - individuals that are crossed out in the raw data
    !(individual %in% c("2525", "2533", "2559", "2560")),
    # Remove any samples that lack a dark measurement (can't measure recovery)
    !(is.na(yield_pre) & is.na(memory_pre)) 
  ) %>%
  # Add column to indicate if sample was a control
  mutate(control = str_detect(notes, regex("control", ignore_case = TRUE))) %>%
  # Fix one entry that had an extra zero
  mutate(
    yield_pre = case_when(
      individual == "2842" & yield_pre > 1000 ~ yield_pre/10,
      TRUE ~ yield_pre
    )) %>%
  mutate(across(contains("memory"), ~ifelse(. < 1, NA, .))) %>%
  assert(within_bounds(0, 1000), contains("yield")) %>%
  assert(within_bounds(1, 4000), contains("memory")) %>%
  mutate(generation = "gametophyte")

# Combine gametophyte data from different years ----
# Also add species names and filter to only filmy ferns

# Read in specimen collection data, which contains species IDs based on
# DNA barcode
# (requires clean_specimen_data.R to be run first)
specimens <- read_csv("data/fern_specimens.csv") 

# Note that when joining collection data to DT data based on individual,
# seven records are dropped because these had collection numbers + subcollection numbers,
# but there is no record of what subcollection number was used in the DT test.
# missing <- bind_rows(gameto_dt_2012, gameto_dt_2013_raw_yield) %>%
#   anti_join(specimens, by = c(individual = "coll_num"))
# 
# specimens %>%
#   filter(str_detect(coll_num, paste(missing$individual, collapse = "|"))) %>%
#   select(genus, specific_epithet, coll_num, date_collected)

gameto_dt_raw_yield <- bind_rows(gameto_dt_2012, gameto_dt_2013_raw_yield) %>% 
  left_join(
    select(specimens, coll_num, species, family),
    by = c(individual = "coll_num")
  ) %>%
  filter(family == "Hymenophyllaceae") %>%
  select(-family) %>%
  mutate(
    generation = "gametophyte",
    control = replace_na(control, FALSE),
    salt = ifelse(control == TRUE, "Control", "MgNO3"),
    dry_time = factor(2, levels = c(2, 15)),
    species = str_replace_all(species, " ", "_")) %>%
  check_dt_data %>%
  # `control` and `file` only apply to gametophyte data
  # - `control`: was this individual used as a control?
  # - `file`: name of file where the raw data come from
  select(
    control, file, species, salt, dry_time, individual, generation, matches("yield|memory|weight")
  )

# Join miniPAM data ----

filmy_dt_raw_yield <-
bind_rows(
  gameto_dt_raw_yield,
  filmy_sporo_dt_raw_yield
)

# Read in data with "fixed" memory values after manually comparing with
# miniPAM data

# This was originally produced by writing out gameto_dt_2013_long: 
# gameto_dt_2013_long %>% 
#   select(individual, file, condition, yield = fvfm, memory, control) %>%
#   write_csv("intermediates/gameto_dt_2013_long.csv")
# and combining with miniPam data as the second sheet to fix errors in 
# the "memory" columns:
# write_csv(minipam_2013_long, "intermediates/minipam_2013_long.csv")
gameto_fixed_lookup <- read_excel("data_raw/intermediates/gameto_dt_2013_long_fix.xlsx", na = c("", "NA", "#VALUE!")) %>%
  select(individual, control, memory, fix, manual_mem, pam_set, note) %>%
  mutate(fix = replace_na(fix, 0)) %>%
  filter(!is.na(memory)) %>%
  mutate(corrected_mem = ifelse(is.na(manual_mem), memory + fix, manual_mem)) %>%
  assert(not_na, corrected_mem) %>%
  select(individual, control, memory, corrected_mem, pam_set, note) %>%
  mutate(generation = "gametophyte")

sporo_fixed_lookup <- read_excel(
  "data_raw/intermediates/sporo_dt_long_fix.xlsx", 
  na = c("", "NA", "#VALUE!")) %>%
  mutate(dry_time = factor(dry_time, levels = c(2, 15))) %>%
  select(species:condition, memory, fix, manual_mem, pam_set, note) %>%
  mutate(fix = replace_na(fix, 0)) %>%
  filter(!is.na(memory)) %>%
  mutate(corrected_mem = ifelse(is.na(manual_mem), memory + fix, manual_mem)) %>%
  assert(not_na, corrected_mem) %>%
  select(species:condition, memory, corrected_mem, pam_set, note) %>%
  mutate(generation = "sporophyte", individual = as.character(individual))

fixed_lookup <- bind_rows(gameto_fixed_lookup, sporo_fixed_lookup)

### Join to miniPAM data to add yield and timestamp ###
# convert to long format
filmy_dt_raw_yield_long <-
filmy_dt_raw_yield %>%
  pivot_longer(names_to = "var_condition", values_to = "value", matches("yield|memory|weight")) %>%
  separate(var_condition, c("var", "condition"))

# Add pam data: gametophytes (do separately from sporophytes because conditions
# for joining are different)
gameto_dt_pam_yield <-
filmy_dt_raw_yield_long %>%
  filter(var == "memory", generation == "gametophyte", !is.na(value)) %>%
  left_join(
    select(fixed_lookup, generation, individual, control, memory, corrected_mem, pam_set, note),
    by = c("generation", "individual", "control", value = "memory")) %>%
  rename(memory = corrected_mem) %>%
  # Join PAM yields based on memory and pam set
  left_join(
    rename(minipam_2013_long, yield_pam = yield), 
    by = c("memory", "pam_set")) %>%
  select(species, salt, dry_time, individual, generation, condition, yield_pam, date_time, note) %>%
  # Make sure no rows got duplicated during the join
  assert_rows(col_concat, is_uniq, individual, salt, condition)

# Add pam data: sporophytes
sporo_dt_pam_yield <-
filmy_dt_raw_yield_long %>%
  filter(var == "memory", generation == "sporophyte", !is.na(value)) %>% # removes 2012 data, which had no memory values recorded
  select(-control, -file) %>%
  left_join(
    select(fixed_lookup, -control),
    by = c("generation", "species", "individual", "salt", "dry_time", "condition", "dataset", value = "memory")) %>%
  rename(memory = corrected_mem) %>%
  # Join PAM yields based on memory and pam set
  left_join(
    rename(minipam_2013_long, yield_pam = yield), 
    by = c("memory", "pam_set")) %>%
  select(species, salt, dry_time, individual, generation, condition, yield_pam, date_time, note, dataset) %>%
  # Make sure no rows got dupliated during the join
  assert_rows(col_concat, is_uniq, species, salt, dry_time, individual, condition)

filmy_dt_pam_yield <- bind_rows(sporo_dt_pam_yield, gameto_dt_pam_yield)

# Select final yields, write out ----
filmy_dt_wide <-
  filmy_dt_raw_yield_long %>%
  pivot_wider(names_from = "var", values_from = "value") %>%
  select(-control, -file) %>%
  rename(yield_manual_entry = yield) %>%
  left_join(
    filmy_dt_pam_yield, 
    by = c("species", "salt", "dry_time", "individual", "generation", "dataset", "condition")) %>%
  mutate(note = replace_na(note, "none")) %>% # so str_detect works
  # Remove rows with no yield or weight at all
  filter(!(is.na(yield_manual_entry) & is.na(yield_pam) & is.na(weight))) %>%
  filter(str_detect(note, "exclude", negate = TRUE)) %>%
  # Convert yield to how it appears on the miniPAM (1000*actual yield)
  mutate(yield_pam = 1000*yield_pam) %>%
  # Check for problems in yield differences between 
  # manual entry and values matched from miniPAM
  mutate(yield_diff = abs(yield_pam - yield_manual_entry)) %>%
  # Flag the check as FALSE if there is a large difference
  # between values entered manually and looked up from miniPAM,
  # and if the note doesn't say to use the miniPAM value.
  mutate(
    pass_check = ifelse(yield_diff > 20 & str_detect(note, "use manual|use pam", negate = TRUE), FALSE, TRUE) %>%
      replace_na(TRUE)) %>%
  assert(isTRUE, pass_check) %>%
  select(-yield_diff, -pass_check) %>%
  # Select the final value to use for yield (manual or miniPAM)
  mutate(yield = case_when(
    # - use manual when indicated in notes
    str_detect(note, "use manual") ~ yield_manual_entry,
    # - use miniPAM when indicated in notes
    str_detect(note, "use pam") ~ yield_pam,
    # - some values only have manual entry, with no memory entered
    !is.na(yield_manual_entry) & is.na(yield_pam) ~ yield_manual_entry,
    # - use the miniPAM otherwise
    TRUE ~ yield_pam
  )) %>%
  select(-yield_manual_entry, -yield_pam, -memory, -note) %>%
  rename(time = date_time) %>%
  # Convert to wide format
  pivot_wider(names_from = "condition", values_from = c("yield", "time", "weight")) %>%
  # Rearrange columns
  select(
    species:dataset,
    yield_pre,
    yield_desiccated,
    yield_30min,
    yield_24hr,
    yield_48hr,
    yield_72hr,
    yield_dry,
    weight_pre,
    weight_desiccated,
    weight_30min,
    weight_24hr,
    weight_48hr,
    weight_72hr,
    weight_dry,
    time_pre,
    time_desiccated,
    time_30min,
    time_24hr,
    time_48hr,
    time_72hr,
    time_dry) %>%
  janitor::remove_empty("cols")

write_csv(filmy_dt_wide, "data/filmy_dt.csv")
  