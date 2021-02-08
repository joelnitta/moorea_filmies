# Clean desiccation chamber raw data (temperature and rel. humidity)

# Load packages
source(here::here("R/packages.R"))

# Load functions
source("R/functions.R")

# 2012 data were measured with __
# FIXME: check name of datalogger
# They lack serial numbers.
# 2012 sporophytes
sporo_2012_chamber <- bind_rows(
  parse_logger_dat("data_raw/2012/filmyDT2_LiCl.csv") %>% mutate(salt = "LiCl"),
  parse_logger_dat("data_raw/2012/filmyDT2_MgNO3.csv") %>% mutate(salt = "MgNO3"),
  parse_logger_dat("data_raw/2012/filmyDT2_NaCl.csv") %>% mutate(salt = "NaCl")
) %>%
  mutate(generation = "sporophyte", year = 2012, species = "mixed") %>%
  # It appears that data before 2012-07-31 18:10:00 UTC was prelminary (NaCl datalogger
  # stabilizes to a different RH than 18%, 58%, or 80%). 
  # Consider the period after this to be the correct data
  filter(date_time > "2012-07-31 18:10:00 UTC")
  
# 2012 gametophytes
gameto_2012_chamber <-
  bind_rows(
    parse_logger_dat("data_raw/2012/gameto1_7-14.csv"),
    parse_logger_dat("data_raw/2012/gameto2_7-24.csv")
  ) %>% mutate(salt = "MgNO3", generation = "gametophyte", year = 2012, species = "mixed")

# 2013 load metadata for hobo loggers
hobo_metadata <- read_csv(
"data_raw/intermediates/hobo_metadata.csv", 
col_types = cols(
  file = col_character(),
  year = col_double(),
  salt = col_character(),
  use_status = col_character(),
  generation = col_character(),
  note = col_character()
))

# Load all data measured with Hobo dataloggers (2013, 2014)
hobo_chamber_data <-
  tibble(
    file = c(
      list.files("data_raw/2013/hobo", pattern = "csv", full.names = TRUE),
      list.files("data_raw/2014/hobo", pattern = "csv", full.names = TRUE)
    )
  ) %>%
  # Add metadata
  left_join(select(hobo_metadata, -serial_no), by = "file") %>%
  mutate(data = map(file, read_hobo)) %>%
  unnest(data) %>%
  # Exclude data not marked "keep" in metadata
  filter(use_status == "keep") %>%
  select(-use_status, -note) %>%
  # modify generation for Hobo s/n 10140659, which was used partly for sporos, partly for gametos
  mutate(
    generation = case_when(
      serial_no == "10140659" & date_time > "2013-07-03" ~ "gametophyte",
      serial_no == "10140659" & date_time < "2013-07-03" ~ "sporophyte",
      TRUE ~ generation
    )
  )

# Combine and write out
combined_chamber_data <- bind_rows(sporo_2012_chamber, gameto_2012_chamber, hobo_chamber_data) %>%
  assert(not_na, date_time, temp, rh, salt, generation, year, species)

write_csv(combined_chamber_data, "data/filmy_dt_chamber.csv")
