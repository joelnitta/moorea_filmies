source("_targets.R")
library(patchwork)

hobo_metadata <- read_csv(
  "data_raw/intermediates/hobo_metadata.csv", 
  col_types = cols(
    file = col_character(),
    year = col_double(),
    species = col_character(),
    duration = col_double(),
    generation = col_character(),
    salt = col_character(),
    use_status = col_character(),
    note = col_character()
  ))

chamber_data <-
tibble(
  file = c(
    list.files("data_raw/2013/hobo", pattern = "csv", full.names = TRUE),
    list.files("data_raw/2014/hobo", pattern = "csv", full.names = TRUE)
  )
) %>%
  left_join(hobo_metadata, by = "file") %>%
  mutate(data = map(file, read_hobo)) %>%
  unnest(data) %>%
  filter(use_status == "keep")

filmy_gameto_dt = load_gameto_dt("data/filmy_gameto_dt.csv")

gameto_times <-
filmy_gameto_dt %>%
  select(individual, contains("time_")) %>%
  pivot_longer(names_to = "condition", values_to = "date_time", -individual) %>%
  filter(!is.na(date_time)) %>%
  mutate(condition = str_remove_all(condition, "time_") %>%
           str_replace_all("30min", "post")) %>%
  filter(condition %in% c("pre", "post")) %>%
  add_count(individual) %>%
  filter(n > 1) %>%
  select(-n) %>%
  pivot_wider(names_from = "condition", values_from = "date_time")

gameto_times_min <- min(gameto_times$pre)
gameto_times_max <- max(gameto_times$post)

dt_times_gameto_2013_plot <-
ggplot(gameto_times) +
  geom_segment(aes(x = pre, xend = post, y = individual, yend = individual))

chamber_gameto_2013 <-
  bind_rows(
chamber_data %>%
  filter(serial_no == "10140659") %>%
  filter(date_time > "2013-07-03"),
chamber_data %>%
  filter(serial_no == "10342564")
)

# Gametophytes 2013
chamber_gameto_2013_rh_plot <-
ggplot(chamber_gameto_2013, aes(x = date_time, y = rh)) +
  geom_line() +
  scale_x_datetime(limits = c(gameto_times_min, gameto_times_max))

chamber_gameto_2013_temp_plot <-
  ggplot(chamber_gameto_2013, aes(x = date_time, y = temp)) +
  geom_line() +
  scale_x_datetime(limits = c(gameto_times_min, gameto_times_max))

chamber_gameto_2013_temp_plot + chamber_gameto_2013_rh_plot + dt_times_gameto_2013_plot + plot_layout(ncol = 1)

# Others...

# SN 10140659 was used for filmies set 1 MgNO3 from 2013-06-16 to 2013-07-02,
# then for gametophytes from 2013-07-05 to 2013-07-21
#
# SN 10342564 was used for gametophytes (despite the name) from 2013-06-09
# to 2013-07-05. The excluded file "data_raw/2013/hobo/Filmies_Mouaputa_200,_400_gameto_MgNO3_2d.csv"
# contains the same set of measurements over a small range of dates (2013-06-09 to 2013-06-11)
chamber_data %>% 
    filter(serial_no %in% c("10342564", "10140659")) %>%
    ggplot(aes(x = date_time, y = rh, color = file)) +
    geom_line() +
    facet_wrap(vars(serial_no), nrow =2, ncol = 1)

# check contents of files
chamber_data %>% 
  filter(date_time < "2014-01-01") %>%
  ggplot(aes(x = date_time, y = rh, color = file, linetype = salt)) +
  geom_line() +
  facet_wrap(vars(serial_no), ncol = 1)