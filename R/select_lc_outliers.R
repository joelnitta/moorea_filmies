# Select outliers in light curve data
source("R/packages.R")

# Load packages for interactive plots/dataframes
library(plotly)
library(crosstalk)
library(DT)
options(persistent = FALSE)

# Define function for interactive plot with dataframe of selected points
# see SO post:
# https://stackoverflow.com/questions/50765687/return-datapoints-selected-in-a-plotly-scatterplot
select_lc_points <- function(plotly_data_tibble, row_select)  {
  
  # Subset data to a single row
  plotly_data <- plotly_data_tibble$data[[row_select]]
  
  plotly_title <- plotly_data_tibble$id[[row_select]]
  
  # Create shared data object so it's in sync between DataTable and Plotly
  shared_data <- crosstalk::SharedData$new(plotly_data)
  
  # Create plot
  plotly_plot <- plotly::plot_ly(shared_data, x = ~par, y = ~etr) %>% 
    plotly::add_lines(y = ~etr_fit) %>%
    plotly::add_markers(alpha = 0.5) %>%
    plotly::highlight("plotly_selected", dynamic = TRUE)
  
  # Create data table
  plotly_data <- DT::datatable(shared_data, extensions = 'Buttons', options = list(
    dom = 'Bfrtip',
    # include download button with file name set to data ID so it can be
    # easily joined later
    buttons = list(
      list(extend = 'csv', filename = plotly_title)),
    text = 'Download'
  ))
  
  # Render the interactive plot
  crosstalk::bscols(widths = c(7, 3), plotly_plot, plotly_data)
  
}

# Load raw light-curve data
filmy_lc <- read_csv("data/filmy_light_curves.csv", col_types = cols(
  type = col_character(),
  no = col_double(),
  f = col_double(),
  fm = col_double(),
  par = col_double(),
  yield = col_double(),
  etr = col_double(),
  date_time = col_datetime(format = ""),
  species = col_character(),
  individual = col_character(),
  generation = col_character(),
  coll_num = col_character(),
  sporo_site = col_character(),
  condition = col_character(),
  date = col_date(format = "")
)
)

# Fit light curves for each individual
filmy_lc_model_fits <-
  filmy_lc %>%
  select(species, generation, individual, coll_num, condition, date, par, etr) %>%
  nest(data = c(par, etr)) %>%
  mutate(
    nls_mod = map(
      data, 
      # nonlinear least-squares estimates of parameters of the subsetted data
      ~nls(
        etr~max(etr)*(1-exp(-k*par)),
        start = list(k = 0.04),
        data =.,
        trace = FALSE,
        control = list(maxiter = 500))
    ),
    k_stats = map(nls_mod, broom::tidy)
  ) %>%
  mutate(fitted = map(nls_mod, ~fitted(.) %>% tibble(etr_fit = .))) %>%
  select(-nls_mod, -k_stats) %>%
  unnest(cols = c(data, fitted)) %>%
  # Convert genus to just first letter
  separate(species, c("genus", "epithet")) %>%
  mutate(genus = substr(genus, 1, 1)) %>%
  unite("species", genus, epithet) %>%
  mutate(id = paste3(species, coll_num, individual, as.character(date), sep = "_")) %>%
  # Nest fitted data
  group_by(species, generation, individual, coll_num, condition, date, id) %>%
  nest() %>%
  ungroup

# Subset to gametophytes
gameto_lc_model_fits <-
  filmy_lc_model_fits %>%
  filter(generation == "gametophyte")

# In RStudio, interactively go through data by increasing row selected one at a time,
# select outliers, and save each as CSV
select_lc_points(gameto_lc_model_fits, 1)
