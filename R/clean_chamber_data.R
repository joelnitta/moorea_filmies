# Clean desiccation chamber raw data (temperature and rel. humidity)

# Load packages
source(here::here("R/packages.R"))

# Load functions
source("R/functions.R")

# Load metadata for data loggers
logger_metadata <- read_csv(
  "data_raw/intermediates/logger_metadata.csv", 
  col_types = cols(
    file = col_character(),
    year = col_double(),
    salt = col_character(),
    use_status = col_character(),
    generation = col_character(),
    note = col_character(),
    serial_no = col_character()
  ))

# 2012 data were measured with Track-It RH/Temp dataloggers for most samples,
# and Hobo ProV2 dataloggers for Callistopteris apiifolia only
# 2012 sporophytes
sporo_2012_chamber <- bind_rows(
  parse_logger_dat("data_raw/2012/track-it/filmyDT2_LiCl.csv") %>% mutate(salt = "LiCl"),
  parse_logger_dat("data_raw/2012/track-it/filmyDT2_MgNO3.csv") %>% mutate(salt = "MgNO3"),
  parse_logger_dat("data_raw/2012/track-it/filmyDT2_NaCl.csv") %>% mutate(salt = "NaCl")
) %>%
  mutate(generation = "sporophyte", year = 2012, species = "mixed") %>%
  # Data before 2012-07-31 18:10:00 UTC was prelminary (NaCl datalogger
  # was actually placed in MgCl2, which wasn't used in final analysis). 
  # Filter to the period after this.
  filter(date_time > "2012-07-31 18:10:00 UTC") %>%
  left_join(
    select(logger_metadata, salt, generation, year, species, file, serial_no),
    c("salt", "generation", "year", "species")
  )
  
# 2012 gametophytes
gameto_2012_chamber <-
  bind_rows(
    parse_logger_dat("data_raw/2012/track-it/gameto1_7-14.csv"),
    parse_logger_dat("data_raw/2012/track-it/gameto2_7-24.csv")
  ) %>% mutate(salt = "MgNO3", generation = "gametophyte", year = 2012, species = "mixed") %>%
  left_join(
    select(logger_metadata, salt, generation, year, species, file, serial_no),
    c("salt", "generation", "year", "species")
  )

# Load all data measured with Hobo dataloggers (2013, 2014)
hobo_chamber_data <-
  tibble(
    file = c(
      list.files("data_raw/2012/hobo", pattern = "csv", full.names = TRUE),
      list.files("data_raw/2013/hobo", pattern = "csv", full.names = TRUE),
      list.files("data_raw/2014/hobo", pattern = "csv", full.names = TRUE)
    )
  ) %>%
  # Add metadata
  left_join(select(logger_metadata, -serial_no), by = "file") %>%
  mutate(data = map(file, read_hobo)) %>%
  unnest(data) %>% 
  # Exclude C. apiifolia NaCl 2012 after 2012-08-02 15:30:00 (end of DT test)
  mutate(
    use_status = case_when(
      species == "Callistopteris_apiifolia" & salt == "NaCl" & year == 2012 & date_time > "2012-08-02 15:30:00" ~ "exclude",
      TRUE ~ use_status
    )
  ) %>%
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
