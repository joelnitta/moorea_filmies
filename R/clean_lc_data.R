# Clean light-curve raw data

# Load packages
library(readxl)
source("R/packages.R")

# Load functions
source("R/functions.R")

# Read in specimen collection data
# (requires clean_specimen_data.R to be run first)
specimens <- read_csv("data/fern_specimens.csv") 

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
  mutate(condition = "lab")

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
field_lc_data <-
field_lc_files %>%
  mutate(lc_data = pmap(., parse_pam)) %>%
  mutate(
    file_short = fs::path_file(file),
    sporo_site = str_replace_all(file_short, "Three Pines 6-30-13", "Three Pines 200m 6-30-13") %>%
      str_match("(^.*m) ") %>% 
      magrittr::extract(,2),
    date = str_match(file_short, " ([0-9]-.*)\\.pam")%>% 
      magrittr::extract(,2) %>%
      mdy) %>%
  select(lc_data, sporo_site, date) %>%
  unnest(cols = c(lc_data)) %>%
  mutate(condition = "field")

# Filter to filmy fern sporophytes
filmy_genera <- c("Abrodictyum", "Callistopteris", "Crepidomanes", "Didymoglossum", "Hymenophyllum", "Polyphlebium", "Vandenboschia")

filmy_sporo_lc_field <-
field_lc_data %>%
  filter(str_detect(id, paste0(filmy_genera, collapse =  "|"))) %>%
  filter(!str_detect(id, regex("fail", ignore_case = TRUE))) %>%
  filter(!str_detect(id, regex("gameto", ignore_case = TRUE))) %>%
  # Exclude Crepidomanes minutum, because didn't record variety and have no way to verify
  filter(!str_detect(id, "minutum")) %>%
  mutate(id = str_remove_all(id, "2422_|_sporo")) %>%
  separate(id, c("genus", "epithet", "individual"), fill = "right", extra = "drop") %>%
  unite("species", genus, epithet) %>%
  mutate(generation = "sporo")
  
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
  select(-family) 

# Combine data and write out ----
filmy_lc_data <- bind_rows(filmy_sporo_lc_lab, filmy_sporo_lc_field, filmy_gameto_lc_field) %>% 
  assert(not_na, type:generation, condition) %>%
  mutate(date = date(date_time))

write_csv(filmy_lc_data, "data/filmy_light_curves.csv")
