# Clean desiccation chamber raw data

# Load packages
source(here::here("R/packages.R"))

# Load functions
source("R/functions.R")

# Parse data from datalogger inside desiccation chamber for sporophytes in 2012
sporo_2012_chamber <- bind_rows(
  parse_logger_dat("data_raw/2012/filmyDT2_LiCl.csv") %>% mutate(salt = "LiCl"),
  parse_logger_dat("data_raw/2012/filmyDT2_MgNO3.csv") %>% mutate(salt = "MgNO3"),
  parse_logger_dat("data_raw/2012/filmyDT2_NaCl.csv") %>% mutate(salt = "NaCl")
)
  
ggplot(sporo_2012_chamber, aes(x = time, y = humidity, color = salt)) +
  geom_line()

ggplot(sporo_2012_chamber, aes(x = time, y = temp, color = salt)) +
  geom_line()

# Parse data from datalogger inside desiccation chamber for gametophytes in 2012
gameto_2012_chamber <-
  bind_rows(
    parse_logger_dat("data_raw/2012/gameto1_7-14.csv"),
    parse_logger_dat("data_raw/2012/gameto2_7-24.csv")
  ) %>% mutate(salt = "MgNO3")

ggplot(gameto_2012_chamber, aes(x = time, y = humidity, color = salt)) +
  geom_line()
