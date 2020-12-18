# Clean specimen collection data

# Load packages
source("R/packages.R")

# Load functions
source("R/functions.R")

# Read in Pteridophyte Phylogeny Group I taxonomic system
# (to add family level taxonomy)
ppgi <- read_csv("data_raw/ppgi_taxonomy_mod.csv")

# Read in collections
fern_specimens_gps_raw <- read_csv("data_raw/specimens.csv") %>%
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
    date_collected = paste(year, month, day, sep = "-")
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
         country, county, locality, site,
         elevation, latitude, longitude,
         observations,
         generation,
         gameto_net_location, gameto_square,
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
