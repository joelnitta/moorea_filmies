# Clean specimen collection data

# Load packages
source(here::here("R/packages.R"))

# Load functions
source("R/functions.R")

# Read in Pteridophyte Phylogeny Group I taxonomic system
# (to add family level taxonomy)
ppgi <- read_csv("data_raw/ppgi_taxonomy_mod.csv", col_types = cols(
  class = col_character(),
  order = col_character(),
  suborder = col_character(),
  family = col_character(),
  subfamily = col_character(),
  genus = col_character(),
  notes = col_character()
))

# Specify column types for specimen data ----
specimens_raw_col_spec <- 
  cols(
    specimenID = col_double(),
    collector_firstname = col_character(),
    collector_lastname = col_character(),
    collector_fullname = col_logical(),
    collectors_other = col_character(),
    collection_number = col_double(),
    subcollection_number = col_character(),
    chromosome_voucher = col_character(),
    coll_no_as_text = col_character(),
    specimen = col_character(),
    genus = col_character(),
    specific_epithet = col_character(),
    infraspecific_rank = col_character(),
    infraspecific_name = col_character(),
    informal_taxon = col_character(),
    morphospecies = col_character(),
    certainty = col_character(),
    country = col_character(),
    county = col_character(),
    province = col_character(),
    locality = col_character(),
    site = col_character(),
    site_full = col_logical(),
    plot = col_character(),
    `Loc Summary` = col_double(),
    observations = col_character(),
    elevation_m = col_double(),
    elevation_error = col_character(),
    elevation_unit = col_double(),
    latitude = col_double(),
    longitude = col_double(),
    date_collected = col_character(),
    date_determined = col_double(),
    determined_by = col_character(),
    created = col_logical(),
    notes = col_character(),
    notes2 = col_character(),
    edits = col_character(),
    modified = col_logical(),
    herbaria = col_character(),
    `Sent To` = col_logical(),
    number_of_sheets = col_character(),
    number_of_duplicates = col_double(),
    `Gameto DT` = col_double(),
    `Gameto Habit` = col_character(),
    `Gameto Net Location` = col_character(),
    `Gameto Square` = col_character(),
    hybrid_formula = col_logical(),
    id_method = col_character(),
    in_plot = col_double(),
    is_gametophyte = col_double(),
    missing_sample = col_double(),
    no_voucher_remaining = col_double(),
    not_a_fern = col_double(),
    `Preserve Type` = col_character(),
    spore_voucher = col_double()
  )


# Read in collections ----
fern_specimens_gps_raw <- read_csv("data_raw/specimens.csv", col_types = specimens_raw_col_spec) %>%
  clean_names() %>%
  # Filter to Moorea or Tahiti collections of ferns by Nitta
  filter(
    (str_detect(collector_lastname, "Nitta") | str_detect(collectors_other, "Nitta")), 
    county %in% c("Moorea", "Tahiti"), 
    is.na(not_a_fern)) %>% 
  # Manipulate columns
  rename(elevation = elevation_m, other_collectors = collectors_other) %>%
  mutate(
    coll_num = paste3(collection_number, subcollection_number, sep = ""),
    specimen = paste3(collector_lastname, coll_num),
    collector = paste(collector_firstname, collector_lastname),
    genus = na_if(genus, "-"),
    specific_epithet = na_if(specific_epithet, "-"),
    species = paste3(genus, specific_epithet),
    taxon = paste3(genus, specific_epithet, infraspecific_name),
    notes = replace_na(notes, ""),
    generation = case_when(
      is_gametophyte == 0 ~ "sporophyte",
      is_gametophyte == 1 ~ "gametophyte",
      TRUE ~ NA_character_
    )
  ) %>%
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
    date_collected = paste3(year, month, day, sep = "-")
  ) %>%
  # Add taxon
  mutate(taxon = paste3(genus, specific_epithet, infraspecific_name)) %>%
  # Add family
  left_join(select(ppgi, genus, family), by = "genus") %>%
  # Select variables
  select(specimen_id, specimen, collector, coll_num,
         family, genus, specific_epithet, infraspecific_rank, infraspecific_name, certainty,
         species, taxon, 
         # scientific_name, author, var_author,
         country, county, locality, site = plot,
         elevation, latitude, longitude,
         observations,
         generation,
         gameto_net_location, gameto_square, gameto_habit,
         date_collected,
         other_collectors,
         herbaria)

# When the data were imported to openrefine, all of the degrees, minutes, and decimal seconds
# got stripped from latitudes and longitudes, so they just appear as normal numbers.
# These are mixed up with decimal degree GPS data.
# Fix them by pulling out those with latitude > 100, which were DMS originally, 
# then converting to decimal degree.
fern_specimens_gps_fix <-
  fern_specimens_gps_raw %>%
  filter(latitude > 100) %>%
  mutate(
    across(c(latitude, longitude), as.character),
    lat_degree = substr(latitude, 1, 2),
    lon_degree = substr(longitude, 1, 3),
    lat_minute = ifelse(
      # Some lat / longs were originally recorded in a *combination* of dms and decimal:
      # e.g., 17°32.881' S 149°50.664'W
      # others normal dms, e.g.  17°32"88' S 149°50"66'W
      # Treat these by detecting period (decimal)
      str_detect(latitude, "\\."), 
      substr(latitude, 3, nchar(latitude)),
      substr(latitude, 3, 4)),
    # In the data for French Polynesia, there are two digits for longitude degree (17 degreess)
    # and three digits for latitude degree (149 degrees)
    lon_minute = ifelse(
      str_detect(longitude, "\\."), 
      substr(longitude, 4, nchar(longitude)),
      substr(longitude, 4, 5)),
    lat_second = ifelse(
      str_detect(latitude, "\\."), 
      NA,
      substr(latitude, 5, 6)),
    lon_second = ifelse(
      str_detect(longitude, "\\."), 
      NA,
      substr(longitude, 6, 7)),
  ) %>%
  mutate(
    across(contains("degree"), ~ifelse(!is.na(.), glue("{.}d"), NA)),
    across(contains("minute"), ~ifelse(!is.na(.), glue("{.}m"), NA)),
    across(contains("second"), ~ifelse(!is.na(.), glue("{.}s"), NA)),
    lat_decimal = paste3(lat_degree, lat_minute, lat_second, "S") %>% sp::char2dms(chm = "m", chs = "s") %>% as.numeric,
    lon_decimal = paste3(lon_degree, lon_minute, lon_second, "W") %>% sp::char2dms(chm = "m", chs = "s") %>% as.numeric
  ) %>%
  select(-c(latitude, longitude, lat_degree, lat_minute, lat_second, lon_degree, lon_minute, lon_second)) %>%
  rename(latitude = lat_decimal, longitude = lon_decimal)
  
# Add back in the specimens with corrected lat/long
fern_specimens <-
  fern_specimens_gps_raw %>%
  anti_join(fern_specimens_gps_fix, by = "specimen_id") %>%
  bind_rows(fern_specimens_gps_fix) %>%
  arrange(specimen_id)

write_csv(fern_specimens, "data/fern_specimens.csv")
