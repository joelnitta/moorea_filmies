# Load packages
source(here::here("R/packages.R"))

# Custom funcs ---
source(here("R/functions.R"))

# Gametophyte 2012 DT data ----

# Read in each excel sheet as a separate dataframe.
# In 2012 only, there were three replicate measurements taken of each individual
# per treatment
gameto_2012_dark <- read_excel("data/2012/all gametos 10-3.xlsx", sheet = "DARK", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x9, yield_3 = x13)

gameto_2012_dry <- read_excel("data/2012/all gametos 10-3.xlsx", sheet = "DRY", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

gameto_2012_30min <- read_excel("data/2012/all gametos 10-3.xlsx", sheet = "10 MIN", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

gameto_2012_24hr <- read_excel("data/2012/all gametos 10-3.xlsx", sheet = "24 HR", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

gameto_2012_48hr <- read_excel("data/2012/all gametos 10-3.xlsx", sheet = "48 HR", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

gameto_2012_72hr <- read_excel("data/2012/all gametos 10-3.xlsx", sheet = "72 HR", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(plot = x1, individual = x2, yield_1 = x5, yield_2 = x10, yield_3 = x15)

# Do the same for an additional set of data from the 200m terrestrial site

gameto_2012_200m_ter_dark <- read_excel("data/2012/new 200 m ter gameto.xlsx", sheet = "Pre", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x7, yield_2 = x12, yield_3 = x17)

gameto_2012_200m_ter_dry <- read_excel("data/2012/new 200 m ter gameto.xlsx", sheet = "Dry + 20 min", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x11, yield_2 = x16, yield_3 = x21)

gameto_2012_200m_ter_30min <- read_excel("data/2012/new 200 m ter gameto.xlsx", sheet = "Dry + 20 min", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x11, yield_2 = x10, yield_3 = x15)

gameto_2012_200m_ter_24hr <- read_excel("data/2012/new 200 m ter gameto.xlsx", sheet = "24 hr", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x7, yield_2 = x12, yield_3 = x17)

gameto_2012_200m_ter_48hr <- read_excel("data/2012/new 200 m ter gameto.xlsx", sheet = "48 hr", skip = 2, col_names = FALSE) %>%
  clean_names() %>%
  select(individual = x1, yield_1 = x7, yield_2 = x12, yield_3 = x17)

gameto_2012_200m_ter_72hr <- read_excel("data/2012/new 200 m ter gameto.xlsx", sheet = "72 hr", skip = 2, col_names = FALSE) %>%
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

# Gametophyte 2013 DT data ----

### Read in gametophyte DT test data ###

# These were originally entered on google drive, and have been
# downloaded from there as xlsx files.
# They have slightly differing column names, but the important
# ones can be identified by containing "yield" or "mem" (memory)

# Read in all gametophyte DT data as a list,
# but don't fix column names yet
gameto_dt_2013_all_raw_names <- 
  # List all 2013 gametophyte DT data sets 2-13 files
  list.files("data/2013", pattern = regex("moorea 2013 gameto"), ignore.case = TRUE, full.names = TRUE) %>%
  # Name by vector of files by the file name, so we can include it when combining the data
  set_names(fs::path_file(.)) %>%
  # Read in the xlsx files as a list
  map(suppressMessages(read_excel)) %>%
  map(clean_names)

# Manually fix some column names that fix_gameto_dt_names() wouldn't handle properly
gameto_dt_2013_all_raw_names[["Moorea 2013 gametophyte DT set1.xlsx"]] <-
  gameto_dt_2013_all_raw_names[["Moorea 2013 gametophyte DT set1.xlsx"]] %>%
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

gameto_dt_2013_all_raw_names[["Moorea 2013 Gameto DT Set 13 DONE.xlsx"]] <-
  gameto_dt_2013_all_raw_names[["Moorea 2013 Gameto DT Set 13 DONE.xlsx"]] %>%
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
gameto_dt_2013_all <- 
  gameto_dt_2013_all_raw_names %>%
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
  assert(within_bounds(1, 4000), contains("memory"))

# Convert to long format
gameto_dt_2013_long <- left_join(
  gameto_dt_2013_all %>%
    select(individual, control, file, contains("yield")) %>%
    pivot_longer(names_to = "condition", values_to = "yield", contains("yield")) %>%
    mutate(condition = str_remove_all(condition, "yield_")),
  gameto_dt_2013_all %>%
    select(individual, control, file, contains("memory")) %>%
    pivot_longer(names_to = "condition", values_to = "memory", contains("memory")) %>%
    mutate(condition = str_remove_all(condition, "memory_")),
  by = c("individual", "control", "file", "condition")
)

# Write out the data in long format for manual comparison with
# miniPAM data to fix errors in the "memory" columns
# gameto_dt_2013_long %>% 
#   select(individual, file, condition, yield = fvfm, memory, control) %>%
#   write_csv("intermediates/gameto_dt_2013_long.csv")
# 
# gameto_dt_2013_long %>%
#   filter(file == "Moorea 2013 Gameto DT Set 13 DONE.xlsx") %>%
#   select(individual, file, condition, yield, memory, control) %>%
#   mutate(condition = str_replace_all(condition, "pre", "dark")) %>%
#   mutate(condition = factor(condition, levels = c("dark", "dry", "30min", "24hr", "48hr", "72hr"), ordered = TRUE)) %>%
#   mutate(file_order = 13) %>%
#   arrange(condition, individual) %>%
#   write_csv("data/intermediates/gameto_dt_2013_long_fix_add_13.csv")


### Read in mini-PAM data

# Loop over the files: read them in, process each data frame in the list, then combine
gameto_fl_2013_long <-
  # List all the 2013 gameto mini pam files
  list.files("data/2013/minipam", pattern = "miniPAM SET\\d ALL.*\\.pam", full.names = TRUE) %>%
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

# Write out the data in long format for manual comparison with
# miniPAM data to fix errors in the "memory" columns
# write_csv(gameto_fl_2013_long, "intermediates/gameto_fl_2013_long.csv")

# Read in data with "fixed" memory values after manually comparing with
# miniPAM data
fixed_lookup <- read_excel("data/intermediates/gameto_dt_2013_long_fix.xlsx", na = c("", "NA", "#VALUE!")) %>%
  select(individual, control, memory, fix, manual_mem, pam_set, note) %>%
  mutate(fix = replace_na(fix, 0)) %>%
  filter(!is.na(memory)) %>%
  mutate(corrected_mem = ifelse(is.na(manual_mem), memory + fix, manual_mem)) %>%
  assert(not_na, corrected_mem) %>%
  select(individual, control, memory, corrected_mem, pam_set, note)

# Map the values from the miniPAM data to the manually entered yield data,
# convert to wide format
gameto_dt_2013 <-
  gameto_dt_2013_long %>%
  # Exclude values with no manually entered yield or memory
  filter(!(is.na(yield) & is.na(memory))) %>%
  rename(yield_manual_entry = yield) %>%
  # Add corrected memory numbers
  left_join(fixed_lookup, by = c("individual", "control", "memory")) %>%
  # DT set 1 is not included in "fixed_lookup" because it was already correct.
  # Fill in 'corrected_mem' for Set 1 as the original memory
  mutate(
    corrected_mem = ifelse(file == "Moorea 2013 gametophyte DT set1.xlsx", memory, corrected_mem),
    pam_set = ifelse(file == "Moorea 2013 gametophyte DT set1.xlsx", "set1", pam_set)) %>%
  select(-memory) %>%
  # Now convert 'corrected_mem' to 'memory'
  rename(memory = corrected_mem) %>%
  # Join PAM yields based on memory and pam set
  left_join(
    rename(gameto_fl_2013_long, yield_pam = yield), 
    by = c("memory", "pam_set")) %>%
  # Convert yield to how it appears on the miniPAM (1000*actual yield)
  mutate(
    yield_pam = 1000*yield_pam,
    note = replace_na(note, "none")) %>%
  filter(str_detect(note, "exclude", negate = TRUE)) %>%
  # Check for problems in yield differences between 
  # manual entry and values matched from miniPAM
  mutate(yield_diff = abs(yield_pam - yield_manual_entry)) %>%
  # flag the check as FALSE if there is a large difference
  # between values entered manually and looked up from miniPAM,
  # and if the note doesn't say to use the miniPAM value.
  mutate(
    pass_check = ifelse(yield_diff > 20 & note != "use pam", FALSE, TRUE) %>%
      replace_na(TRUE)) %>%
  assert(isTRUE, pass_check) %>%
  select(-yield_error, -yield_diff, -pass_check) %>%
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
  select(-yield_manual_entry, -yield_pam) %>%
  # Check all rows have a yield entry
  assert(not_na, yield) %>%
  select(-file, -memory, -pam_set, -note) %>%
  rename(time = date_time) %>%
  # Convert to wide format
  pivot_wider(names_from = "condition", values_from = c("yield", "time")) %>%
  # Add column for desiccation treatment and dry time
  mutate(
    salt = ifelse(isTRUE(control), "Control", "MgNO3"),
    dry_time = 2) %>%
  select(-control)
    
### Combine all datasets ----
gameto_dt <- bind_rows(gameto_dt_2012, gameto_dt_2013)

# read in ppgi
ppgi <- read_csv("data/ppgi_taxonomy_mod.csv")

#### Filter to only filmy ferns ----

# Read in collections
collections <- read_csv("data/specimens.csv") %>%
  clean_names() %>%
  filter(collector_lastname == "Nitta") %>%
  # Manipulate columns
  mutate(coll_num = paste3(collection_number, subcollection_number, sep = "")) %>%
  mutate(specimen = paste3(collector_lastname, coll_num)) %>%
  mutate(collector = paste(collector_firstname, collector_lastname)) %>%
  rename(elevation = elevation_m, other_collectors = collectors_other) %>%
  mutate(species = paste3(genus, specific_epithet)) %>%
  mutate(taxon = paste3(genus, specific_epithet, infraspecific_name)) %>%
  mutate(notes = replace_na(notes, "")) %>%
  # Format collection date YYYY-MM-DD
  mutate(
    date_collected = date_collected %>%
      str_remove_all("^[:alpha:]+ ") %>% 
      str_remove_all("00:00:00 [:alpha:]+ ") %>%
      mdy() %>%
      as_date(),
    year = year(date_collected),
    month = month(date_collected) %>% str_pad(side = "left", pad = "0", width = 2),
    day = day(date_collected) %>% str_pad(side = "left", pad = "0", width = 2),
    date_collected = paste(year, month, day, sep = "-")
  ) %>%
  # Select variables
  select(specimen_id, specimen, coll_num,
         genus, specific_epithet, infraspecific_rank, infraspecific_name, certainty,
         # species, taxon, scientific_name, author, var_author,
         country, locality, site, observations,
         elevation, latitude, longitude,
         collector, other_collectors,
         herbaria,
         date_collected) %>%
  mutate(taxon = paste3(genus, specific_epithet, infraspecific_name)) %>%
  left_join(ppgi, by = "genus")

# Note that when joining collection data to DT data based on individual, some
# records are dropped because these had collection numbers + subcollection numbers,
# but there is no record of what subcollection number was used in th DT test.
missing <- gameto_dt %>%
  anti_join(collections, by = c(individual = "coll_num"))

read_csv("data/specimens.csv") %>%
  clean_names() %>%
  filter(collection_number %in% missing$individual) %>%
  select(genus, specific_epithet, collection_number, subcollection_number, date_collected)

# Add taxon and family, filter to only filmy ferns
filmy_gameto_dt <-
gameto_dt %>%
  left_join(
    select(collections, coll_num, taxon, family),
    by = c(individual = "coll_num")
  ) %>%
  filter(family == "Hymenophyllaceae") %>%
  select(-family)

write_csv(filmy_gameto_dt, "filmy_gameto_dt.csv")
