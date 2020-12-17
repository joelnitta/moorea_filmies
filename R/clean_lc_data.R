# Clean light-curve raw data

# Load packages
library(readxl)
source("R/packages.R")

# Load functions
source("R/functions.R")

# Filmy fern sporophytes ----
filmy_lc_files <- c(
  "data_raw/2013/minipam/Filmy Fern  light curves 8-1-13.pam",
  "data_raw/2013/minipam/Filmy Fern  light curves 8-6-13.pam",
  "data_raw/2013/minipam/Filmy Fern  light curves 8-8-13.pam",
  "data_raw/2013/minipam/Abrodictyum_dentatum light curves 7-31-13.pam"
)

filmy_lc_data <- map_df(filmy_lc_files, ~parse_pam(., ret_type = "lc")) %>% 
  # Filter out non-labeled data, standardize names
  mutate(num_only = str_detect(id, "[a-z]|[A-Z]", negate = TRUE)) %>%
  filter(!is.na(id), !num_only) %>%
  select(-num_only) %>% 
  mutate(
    id = str_remove_all(id, "3090_|3171_| l-int3|\\.") %>%
      str_replace_all("Abrodoctyum_cfasaegrayi", "Abrodictyum_asaegrayi") %>% # FIXME: Need to check which one
      str_replace_all("C_apiifolia", "Callistopteris_apiifolia") %>%
      str_replace_all("H_multifidum", "Hymenophyllum_multifidum") %>%
      str_replace_all("H_palladum", "Hymenophyllum_pallidum") %>%
      str_replace_all("Hymenophyllum_sp", "Hymenophyllum_braithwaitei") %>%
      str_replace_all("V_maxima", "Vandenboschia_maxima")
    ) %>%
  separate(id, c("genus", "epithet", "individual"), remove = FALSE) %>%
  unite("species", genus:epithet) %>%
  select(-id)

# Gametophytes ----

# Some of the rows in the .pam files have errors (e.g., light curve started
# but didn't finish), so manually inspect and exclude these
# so that the data can be read in properly.
gameto_lc_files <- tribble(
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

gameto_lc_data <- pmap_df(gameto_lc_files, parse_pam)
