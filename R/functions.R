# Define functions used in the analysis

# Data preprocessing ----

# This functions are used for pre-processing raw data (scripts `R/clean_*.R`)

#' Read in data from a hobo datalogger.
#' 
#' Hobo dataloggers (Onset Corp.) store data in properietary ".hobo" format.
#' This needs to be convert to CSV first by opening the .hobo file with
#' HOBOware software (https://www.onsetcomp.com/hoboware-free-download/),
#' then exporting as CSV.
#'
#' @param file Path to CSV file exported from HOBOware
#'
#' @return Dataframe including columns:
#' - temperature: temperature (degrees celsius)
#' - rh: relative humidity (percent)
#' - serial_no: serial number
#' 
read_hobo <- function(file) {
  
  # Get original column names
  # these contain the serial number, so will be different for each file
  col_names <-
    suppressMessages(suppressWarnings(readr::read_csv(file))) %>%
    janitor::clean_names() %>%
    colnames()
  
  # Extract serial number 
  serial_no <- col_names %>%
    magrittr::extract(3) %>%
    stringr::str_match("s_n_([0-9]+)") %>%
    magrittr::extract(,2)
  
  # Make a vector of clean column names that can be used
  # for any hobo csv file
  clean_col_names <-
    col_names %>%
    str_match("date|time|temp|rh|coupler_detached|coupler_attached|stopped|host_connected|end_of_file") %>%
    magrittr::extract(,1)
  
  # Read in the csv file using the clean column names
  readr::read_csv(
    file,
    skip = 1,
    col_names = clean_col_names,
    col_types = cols(
      date = col_character(),
      time = col_time(format = ""),
      temp = col_double(),
      rh = col_double(),
      .default = col_character()
    )) %>%
    # parse date_time
    dplyr::mutate(date_time = paste(date, time) %>% lubridate::ymd_hms()) %>%
    # only include date_time, temperature, rel. humidity, and serial number
    dplyr::select(date_time, temp, rh) %>%
    ggplot2::remove_missing(na.rm = TRUE) %>%
    dplyr::mutate(serial_no = serial_no)
  
}

#' Convert a factor to numbers
#'
#' @param factor 
#'
#' @return factor, with all levels replaced with numbers
#' @export
#'
#' @examples
#' numberify(
#' factor(c("a", "a", "b"), levels = c("a", "b", "c"))
#' )
numberify <- function (factor) {
  assertthat::assert_that(is.factor(factor))
  levels(factor) <- 1:length(levels(factor))
  factor
}


#' Verify that filmy fern desiccation tolerance test data are valid
#'
#' @param data Dataframe; filmy fern desiccation tolerance test data 
#'
#' @return dataframe
#' 
check_dt_data <- function(data) {
  
  filmy_species <- c(
    "Abrodictyum_asaegrayi1",
    "Abrodictyum_asaegrayi2",
    "Abrodictyum_caudatum",
    "Abrodictyum_dentatum",
    "Callistopteris_apiifolia",
    "Crepidomanes_bipunctatum",
    "Crepidomanes_humile",
    "Crepidomanes_kurzii",
    "Crepidomanes_minutum1",
    "Crepidomanes_minutum2",
    "Crepidomanes_minutum3",
    "Didymoglossum_tahitense",
    "Hymenophyllum_braithwaitei",
    "Hymenophyllum_digitatum",
    "Hymenophyllum_flabellatum",
    "Hymenophyllum_javanicum",
    "Hymenophyllum_multifidum",
    "Hymenophyllum_pallidum",
    "Hymenophyllum_polyanthos",
    "Polyphlebium_borbonicum",
    "Polyphlebium_endlicherianum",
    "Vandenboschia_maxima")
  
  filmy_salts <- c("Control", "H2O", "LiCl", "MgNO3", "NaCl")
  
  data %>%
    assertr::assert(assertr::in_set(filmy_salts, allow.na = FALSE), salt) %>%
    assertr::assert(assertr::in_set(filmy_species, allow.na = FALSE), species) %>%
    assertr::assert(assertr::in_set(2,15), dry_time) %>%
    assertr::assert(is.factor, dry_time) %>%
    assertr::assert(is.numeric, matches("yield|weight"))
  
}

#' Parse a *.pam file
#' 
#' This is file format output by the miniPAM photosynthesis yield analyzer
#'
#' @param file Path to *.pam file
#' @param exclude_lines Lines to exclude from the raw *.pam file (useful if some
#' of the data got corrupted)
#' @param ret_type Either "lc" for light-curves, or "fl" for fluorescence
#' @param recalc_yield Boolean; should the yield (Fv/Fm) be recalculated from 
#' Fm and F0? (miniPAM returns "-" for the yield if Fm and F0 are low, but we
#' can recalculate it ourselves). In this case, a boolean "yield_error" column will
#' be added to the results to indicate if the original yield calculation would have
#' returned "-".
#'
#' @return Dataframe
parse_pam <- function (file, exclude_lines = NULL, ret_type = "lc", recalc_yield = FALSE) {
  
  assertthat::assert_that(is.character(ret_type))
  
  assertthat::assert_that(ret_type %in% c("lc", "fl"),
                          msg = "'ret_type' must be either 'lc' or 'fl'")
  
  # The .pam data format is basically tab-delimited text (the delimiter is a semicolon),
  # but the number and type of column depends on different data types. Each row has a 
  # single data type. They are as follows:
  #
  # V: wincontrol version
  # C: configuration info
  # D: device info
  # F: Fluorescence measurement
  # FO: Fluorescence measurement, at start of light curve
  # REG1: data (or warning) at start of light curve
  # REG2: data (or warning) at start of light curve
  # SCHS: Chart start
  # SLCE: Light curve end
  # SLCS: Light curve start
  
  # First read in PAM file as raw character vector
  data_pam <- read_lines(file)
  
  # Exclude lines if applicable
  if(!is.null(exclude_lines)) data_pam <- data_pam[-exclude_lines]
  
  # Extract data types by line
  data_types <- map_chr(data_pam, ~str_split(., ";") %>% map_chr(3) %>% str_trim)
  
  # Fluorescence data columns (data types "F", "FO")
  fl_cols <- c("date", "time", "type", "no", "x1", "x2", "f", "fm", "par", "yield", "etr", "x3", "x4", "x5")
  
  # read_delim will issue some warnings because some rows have extra columns including
  # (to me) unknown data that can be ignored
  fl_data <- suppressWarnings(read_delim(data_pam[data_types %in% c("F", "FO")], delim = ";", col_names = fl_cols)) %>%
    mutate(date_time  = lubridate::ymd(date) + lubridate::hms(time)) %>%
    select(-date, -time)
  
  # The miniPAM automatically sets yield to "-" if there was not enough signal in
  # fm and f0. Optionally over-ride this and calculate fv/fm ourselves.
  # In this case, add a flag "yield_error" set to TRUE if the miniPAM measurement
  # would have been an error
  if(isTRUE(recalc_yield)) fl_data <- fl_data %>%
    mutate(yield_recalc = (fm - f)/fm) %>%
    mutate(yield_error = ifelse(yield == "-", TRUE, FALSE)) %>%
    # Treat negative yields as zero
    mutate(yield_recalc = ifelse(yield_recalc < 0, 0, yield_recalc)) %>%
    mutate(yield = yield_recalc) %>%
    select(-yield_recalc)
  
  # Exit early returning only fluorescence data if return type is "fl"
  if(ret_type == "fl") return(select(fl_data, type, memory = no, f, fm, par, contains("yield"), etr, date_time)) 
  
  # Light curve data columns (data types "SLCS", "SLCE")
  # - comment (usually) includes the specimen name, but is sometimes missing (NA)
  lc_cols <- c("date", "time", "type", "status", "comment")
  
  # read_delim will issue some warnings because some rows have extra columns including
  # (to me) unknown data that can be ignored
  
  # goal is to reform LC data so we have one row per specimen, with start and end times.
  lc_data_intermed <- suppressWarnings(read_delim(data_pam[data_types %in% c("SLCS", "SLCE")], delim = ";", col_names = lc_cols)) %>%
    arrange(date, time) %>%
    # verify that LC data are a series of starts/ends. Each pair corresponds to measurement of a specimen.
    verify(.$status == rep(c("Light Curve start", "Light Curve end"), nrow(.)/2)) %>%
    mutate(date_time  = lubridate::ymd(date) + lubridate::hms(time)) %>%
    mutate(status = str_remove_all(status, "Light Curve ")) %>%
    select(status, comment, date_time) %>%
    mutate(temp_id = rep(1:(nrow(.)/2), 2) %>% sort)
  
  specimens <- lc_data_intermed %>%
    filter(!is.na(comment)) %>%
    select(specimen = comment, temp_id)
  
  lc_data <-
    lc_data_intermed %>%
    left_join(specimens, by = "temp_id") %>%
    mutate(id = ifelse(is.na(specimen), temp_id, specimen)) %>%
    pivot_wider(id_cols = id, names_from = "status", values_from = "date_time")
  
  # Join FL and LC data: now we have the specimen linked to the fl data for each light curve
  suppressWarnings(
    fuzzyjoin::fuzzy_left_join(
      fl_data, lc_data,
      by = c(
        "date_time" = "start",
        "date_time" = "end"
      ),
      match_fun = list(`>=`, `<=`)
    ) %>%
      select(type, no, f, fm, par, yield, etr, date_time, id) %>%
      mutate(
        across(where(is.character) & matches("f|fm|par|yield|etr"), parse_number),
        across(c(type, id), as.character))
  )
}

#' Parse datalogger data for a desiccation chamber
#'
#' FIXME: find the name of this type of datalogger
#'
#' @param file Path to datalogger data in CSV format
#'
#' @return Dataframe
parse_logger_dat <- function (file) {
  suppressWarnings(suppressMessages(readr::read_csv(file))) %>%
    janitor::clean_names() %>%
    janitor::remove_empty(which = "cols") %>%
    transmute(
      date_time = mdy_hms(time_mm_dd_yyyy_hh_mm_ss),
      temp = chan_1_deg_c,
      rh = chan_2_percent_rh)
}

#' Standardize column names in gametophyte DT data
#'
#'  @param data A dataframe
#'
#' @return Dataframe
fix_gameto_dt_names <- function (data) {
  
  data %>%
    # There are multiple "notes" columns. Combine these into a single column.
    # "notes" sometimes indicate if a sample was a control!
    rowwise(.) %>%
    mutate(notes_combined = jntools::paste3(c_across(matches("note")), collapse = " ")) %>%
    ungroup() %>%
    # Standardize column names across files.
    # "gametophytes" contains the collection number
    rename_with(., ~"individual", matches("^gametophytes$")) %>%
    # For some reason, the initial measurement (pre-desiccation) was called "light curve"
    # e.g., "Light Curve Yield", "Light Curve Memory"
    rename_with(., ~"memory_pre", matches("light_curve") & matches("_mem")) %>%
    rename_with(., ~"yield_pre", matches("light_curve") & matches("yield")) %>%
    # - memory and yield in the dry condition
    rename_with(., ~"memory_dry", matches("dry") & matches("_mem")) %>%
    rename_with(., ~"yield_dry", matches("dry") & matches("yield")) %>%
    # - memory and yield 30 minutes after recovery
    rename_with(., ~"memory_30min", matches("rewet") & matches("_mem")) %>%
    rename_with(., ~"yield_30min", matches("rewet") & matches("yield")) %>%
    # - memory and yield 24 hr after recovery
    rename_with(., ~"memory_24hr", matches("24_hr") & matches("_mem")) %>%
    rename_with(., ~"yield_24hr", matches("24_hr") & matches("yield")) %>%
    # - memory and yield 48 hr after recovery
    rename_with(., ~"memory_48hr", matches("48_hr") & matches("_mem")) %>%
    rename_with(., ~"yield_48hr", matches("48_hr") & matches("yield")) %>%
    # - memory and yield 72 hr after recovery
    rename_with(., ~"memory_72hr", matches("72_hr") & matches("_mem")) %>%
    rename_with(., ~"yield_72hr", matches("72_hr") & matches("yield")) %>%
    select(., individual, notes = notes_combined, matches("^memory_|^yield_"))%>%
    mutate(., individual = as.character(individual))%>%
    # Force yield and memory columns to numeric
    mutate(., across(-c(individual, notes) & where(is.character), parse_number))
  
}

#' Interactively select points in a scatterplot and save them to a CSV file
#' 
#' For selecting and excluding outliers from light-curve plots.
#' 
#' Based on this SO post:
#' https://stackoverflow.com/questions/50765687/return-datapoints-selected-in-a-plotly-scatterplot
#'
#' @param plotly_data_tibble Nested data frame with one row per species,
#' including an `id` column with the sample ID, and a `data` list-column with
#' nested columns `etr`, `etr_fit`, and `par`
#' @param row_select Index of the nested data frame row that should be plotted
#'
#' @return Nothing; in an RStudio session, the interactive plot will show up in 
#' the Viewer pane.
#' 
select_lc_points <- function(plotly_data_tibble, row_select)  {
  
  # Subset data to a single row
  plotly_data <- plotly_data_tibble %>%
    magrittr::extract(row_select,) %>%
    select(species, generation, individual, coll_num, condition, date, data, fitted) %>%
    unnest(cols = c(data, fitted)) %>%
    select(etr, par, etr_fit, everything())
  
  plotly_title <- plotly_data %>%
    abbrev_sp %>%
    mutate(label = paste(species, generation, individual, coll_num, condition, date, sep = "_")) %>%
    pull(label) %>%
    unique()
  
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

# Data loading ----

#' Load data from a desiccation tolerance (DT) experiment
#'
#' @param file Path to data file (CSV)
#'
#' @return Dataframe
#' 
load_filmy_dt <- function (file) {
  readr::read_csv(file, col_types = cols(
    species = col_character(),
    salt = col_factor(),
    dry_time = col_factor(),
    individual = col_character(),
    dataset = col_character(),
    generation = col_character(),
    yield_pre = col_double(),
    yield_30min = col_double(),
    yield_24hr = col_double(),
    yield_48hr = col_double(),
    yield_72hr = col_double(),
    yield_dry = col_double(),
    weight_pre = col_double(),
    weight_desiccated = col_double(),
    weight_30min = col_double(),
    weight_24hr = col_double(),
    weight_48hr = col_double(),
    weight_dry = col_double(),
    time_pre = col_datetime(format = ""),
    time_30min = col_datetime(format = ""),
    time_24hr = col_datetime(format = ""),
    time_48hr = col_datetime(format = ""),
    time_72hr = col_datetime(format = ""),
    time_dry = col_datetime(format = "")
  )) %>%
    mutate(
      section = case_when(
        species == "Callistopteris_apiifolia" & dataset == "2012" ~ "si",
        species == "Callistopteris_apiifolia" & salt == "H2O" ~ "si",
        species == "Hymenophyllum_polyanthos" & dataset == "2013_4" ~ "si",
        species == "Polyphlebium_borbonicum" & dataset == "2013_2" ~ "si",
        salt == "Control" ~ "si",
        TRUE ~ "main"
      )
    ) %>%
    check_dt_data
}

load_filmy_dt_chamber <- function (file) {
  readr::read_csv(file,
  col_types = cols(
    date_time = col_datetime(format = ""),
    temp = col_double(),
    rh = col_double(),
    salt = col_character(),
    generation = col_character(),
    year = col_double(),
    serial_no = col_character()
  )
  )
}

load_filmy_lc <- function (file) {
  readr::read_csv(file, col_types = cols(
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
}

#' Unzip Nitta et al 2017 Ecol Mono data file downloaded
#' from Dryad (https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g)
#' and extract needed data files.
#' 
#' (This is used to download the file as part of the function, but
#' as it is not clear how to do this with the new Dryad API that
#' was implemented ca. Sept 2019, it now requires the user to have 
#' already downloaded the zipped data file from Dryad manually.)
#' 
#' @param zipped_path Name of downloaded zip file.
#' @param unzip_path Path to directory to put the unzipped
#' contents (will be created if needed).
#' @param ... Extra arguments; not used by this function, but
#' meant for tracking with drake.
#' @return Four unzipped data files:
#' - all_plots.csv: Community matrix for ferns of Moorea and Tahiti
#' including sporophytes (plot names with "_S") and gametophytes
#' (plot names with "_G").
#' - treepl_Moorea_Tahiti.tre: Dated tree for pteridophytes of Moorea and Tahiti
#' - sites.csv: Site metadata including name, latitude, longitude, and elevation (m)
#' - species.csv: List of species included in the Nitta et al. 2017 study
#'
unzip_nitta_2017 <- function (zipped_path, unzip_path, ...) {
  
  # The dryad file is actually a zipped folder within a zipped folder.
  # Unzip the first one to a temporary folder, then unzip the rest from there.
  temp_dir <- tempdir()
  
  unzip(zipped_path, exdir = temp_dir, overwrite = TRUE)
  
  temp_zip <- fs::path(temp_dir, "data_and_scripts.zip")
  
  # Unzip only needed data files to data/nitta_2017/
  unzip(temp_zip, "data_and_scripts/Comm_Phylo_Analysis/data/all_plots.csv", exdir = unzip_path, junkpaths = TRUE, overwrite = TRUE)
  unzip(temp_zip, "data_and_scripts/Comm_Phylo_Analysis/data/treepl_Moorea_Tahiti.tre", exdir = unzip_path, junkpaths = TRUE, overwrite = TRUE)
  unzip(temp_zip, "data_and_scripts/shared_data/sites.csv", exdir = unzip_path, junkpaths = TRUE, overwrite = TRUE)
  
  # Cleanup
  fs::file_delete(temp_zip)
  
}

# Data wrangling ----

#' Cluster rows of a dataframe based on a single column
#' 
#' For example, start times in a dataframe of desiccation tolerance
#' start and end times.
#'
#' @param data Data to cluster
#' @param col Column to use for clustering
#' @param k Number of clusters to generate
#'
#' @return Dataframe with column "cluster" added
#' 
add_clusters <- function(data, col, k) {
  kmeans_uni_groups <- Ckmeans.1d.dp::Ckmeans.1d.dp(as.numeric(data[[col]]), k = k)
  mutate(data, cluster = kmeans_uni_groups$cluster)
}

#' Calculate % recovery from desiccation in a desiccation tolerance (DT) experiment
#'
#' @param data Dataframe; data read in from DT experiment
#'
#' @return Dataframe
calc_recovery <- function (data) {
  data %>%
    # Make sure pre-treatment yield must be greater than zero, or recovery will be Inf
    verify(yield_pre > 0) %>%
    assert(not_na, yield_pre) %>%
    select(species, salt, dry_time, individual, generation, contains("yield")) %>%
    select(-yield_72hr, -contains("yield_dry")) %>%
    pivot_longer(names_to = "rec_time", values_to = "yield_recover", matches("30|24|48")) %>%
    mutate(
      rec_time = str_remove_all(rec_time, "yield_"),
      recovery = yield_recover / yield_pre) %>%
    select(-contains("yield"))
}

#' Calculate mean % recovery from desiccation in a desiccation tolerance (DT) experiment
#'
#' @param data Dataframe; DT data
#'
#' @return Dataframe
calc_mean_recovery <- function (data) {
  data %>%
    mutate(rec_time = factor(rec_time, levels = c("30min", "24hr", "48hr"))) %>%
    group_by(species, salt, dry_time, rec_time, generation) %>%
    summarize(
      recovery_mean = mean(recovery, na.rm = TRUE),
      recovery_sd = sd(recovery, na.rm = TRUE),
      recovery_n = n(),
      .groups = "drop")
}

#' Process raw community matrix 
#'
#' @param community_matrix_path Path to community matrix of Nitta et al 2017
#' @param species_list Vector of species to include
#' @param moorea_sites Dataframe of sites on Moorea in this study
#'
#' @return Tibble (community data matrix)
#' 
process_community_matrix <- function (community_matrix_path, species_list, moorea_sites) {
  read_csv(community_matrix_path) %>%
    rename(site = X1) %>% 
    gather(species, abundance, -site) %>%
    filter(str_detect(site, "_S")) %>%
    mutate(site = str_remove_all(site, "_S")) %>%
    # Keep only sites on Moorea
    filter(site %in% moorea_sites$site) %>%
    # Keep only species in species list (ferns of Moorea)
    filter(species %in% species_list) %>%
    spread(site, abundance)
}



calculate_indiv_etr <- function (data) {
  
  # Calculate maximum ETR for each individual
  data %>%
    assert(not_na, etr) %>%
    group_by(species, generation, individual) %>%
    summarize(
      etr_max = max(etr),
      .groups = "drop"
    )
  
}

calculate_mean_etr <- function (data) {
  
  # Calculate mean ETR max by species and generation
  data %>%
    assert(not_na, etr_max) %>%
    group_by(generation, species) %>%
    summarize(
      etr_mean = mean(etr_max),
      etr_sd = sd(etr_max),
      etr_n = n()
    ) %>%
    ungroup()
  
}

calculate_indiv_par <- function (data) {
  
  # Calculate critical PAR value for each individual
  data %>%
    select(species, generation, individual, par, etr) %>%
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
    select(species, generation, individual, k_stats) %>%
    unnest(cols = c(k_stats)) %>%
    # calculate critical par (par where reach 95% of max etr)
    mutate(par_critical = -log(0.05)/estimate)
  
}

calculate_mean_par <- function (data) {
  
  # Next calculate the mean value by species and generation
  data %>%
    group_by(generation, species) %>%
    summarize(
      par_mean = mean(par_critical),
      par_sd = sd(par_critical),
      par_n = n()
    ) %>%
    ungroup()
  
}

combine_mean_phys_traits <- function (recovery_species_means, etr_species_means, par_species_means) {
  
  recover_filtered <-
    recovery_species_means %>%
    mutate(
      keep = case_when(
        # For sporophytes, use only recovery after 48 hr for 2-day drying treatment with MgNO3
        # But make exception for A. dentatum, which was NaCl
        species == "Abrodictyum_dentatum" & generation == "sporophyte" & salt == "NaCl" & dry_time == "2" & rec_time == "48hr" ~ TRUE,
        salt == "MgNO3" & dry_time == "2" & rec_time == "48hr" ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(keep) %>%
    select(species, generation, starts_with("recover"))
  
  list(
    recover_filtered, etr_species_means, par_species_means
  ) %>%
    reduce(left_join, by = c("species", "generation"))
  
}


#' Calculate mean VPD from climate data
#'
#' @param climate Dataframe with climate data including
#' relative humidity, temperature, and VPD calculated from these measured
#' once every 15 minutes by site and growth habit (terrestrial or epiphytic)
#'
#' @return Tibble with mean VPD by site and growth habit 
calculate_mean_vpd <- function (climate) {
  climate %>% 
    group_by(site, habit) %>% 
    summarize(vpd = mean(vpd), .groups = "drop") %>% 
    assert(not_na, vpd)
}

#' Calculate mean VPD for filmy fern gametophytes
#'
#' Treats specimens collected from rock walls ("epipetric") as epiphytic
#'
#' @param mean_vpd Mean Vapor Pressure Deficit (VPD) measured with dataloggers
#' at each site
#' @param specimens_raw Raw specimen collection data
#' @param moorea_sites Dataframe of collection sites on Moorea
#' @param filmy_species Vector of filmy fern species names
#'
#' @return Tibble with mean VPD per species observed for gametophytes.
#' 
calculate_mean_vpd_gameto <- function(
  mean_vpd,
  specimens_raw,
  moorea_sites,
  filmy_species
) {
  
  gameto_filmy_specimens <-
    specimens_raw %>%
    # only use needed columns
    select(starts_with("gameto"), specimen, genus, specific_epithet, generation, site) %>%
    mutate(
      species = paste(genus, specific_epithet, sep = "_"),
      # **Important**: treat epipteric as epiphytic
      gameto_habit = str_to_lower(gameto_habit) %>% str_replace_all("epipetric", "epiphytic")) %>%
    # subset to filmy fern gametophytes in plots on Moorea
    filter(
      species %in% filmy_species, 
      generation == "gametophyte",
      !is.na(gameto_habit),
      site %in% moorea_sites$site) %>%
    select(specimen, species, site, habit = gameto_habit)
  
  # Calculate mean VPD based on occurrences of gametophytes, including
  # growth habit (treating epipetric as epiphytic)
  gameto_filmy_specimens %>%
    left_join(mean_vpd, by = c("site", "habit")) %>%
    group_by(species) %>%
    # there are several sites that are missing a datalogger,
    # so don't include na values when calculating mean
    summarize(
      gameto_vpd = mean(vpd, na.rm = TRUE),
      .groups = "drop"
    )
  
}

#' Calculate mean VPD for filmy fern sporophytes
#'
#' Treats specimens collected on rocks on the ground surface ("saxicolous") as terrestrial
#'
#' @param community_matrix_raw Raw community matrix data
#' @param mean_vpd Mean Vapor Pressure Deficit (VPD) measured with dataloggers
#' at each site
#' @param traits Growth habit for each species
#' @param moorea_sites Dataframe of collection sites on Moorea
#' @param filmy_species Vector of filmy fern species names
#'
#' @return Tibble with mean VPD per species observed for sporophytes
calculate_mean_vpd_sporo <- function (community_matrix_raw, filmy_species, moorea_sites, traits, mean_vpd) {
  # Convert community matrix to long form
  community_matrix_raw %>%
    rename(site = X1) %>% 
    gather(species, abundance, -site) %>%
    filter(species %in% filmy_species) %>%
    mutate(generation = case_when(
      str_detect(site, "_S") ~ "sporophyte",
      str_detect(site, "_G") ~ "gametophyte"
    )) %>%
    mutate(site = str_remove_all(site, "_S") %>% str_remove_all("_G")) %>%
    # Keep only filmy fern sporophytes in sites on Moorea
    filter(site %in% moorea_sites$site, generation == "sporophyte", abundance > 0) %>%
    select(-generation) %>%
    left_join(select(traits, species, habit), by = "species") %>%
    assert(not_na, habit) %>%
    # For sporophytes, treat "saxicolous" as "terrestrial"
    mutate(habit = case_when(
      str_detect(habit, regex("terr|saxi", ignore_case = TRUE)) ~ "terrestrial",
      str_detect(habit, regex("epi", ignore_case = TRUE)) ~ "epiphytic",
      TRUE ~ NA_character_
    )) %>%
    assert(not_na, habit) %>%
    left_join(mean_vpd, by = c("site", "habit")) %>%
    filter(!is.na(vpd)) %>%
    # To get mean VPD by species x habit, "uncount" abundances to one occurrence per row
    # (like gametophytes)
    uncount(abundance) %>%
    group_by(species) %>%
    summarize(
      sporo_vpd = mean(vpd),
      .groups = "drop"
    )
}

#' Calculate relative water content for sporophytes at the individual level
#'
#' @param data Filmy fern desiccation tolerance test data
#'
#' @return Dataframe
#' 
calculate_indiv_water <- function (data) {
  data %>%
    filter(generation == "sporophyte") %>%
    select(-generation) %>% 
    select(species:dataset, contains("weight")) %>%
    pivot_longer(names_to = "rec_time", values_to = "weight", weight_desiccated:weight_48hr) %>%
    mutate(rel_water_content =  (weight - weight_dry) / (weight_pre - weight_dry)) %>%
    filter(!is.na(rel_water_content)) %>%
    mutate(
      rec_time = str_remove_all(rec_time, "weight_") %>%
        factor(levels = c("desiccated", "30min", "24hr", "48hr"))) %>%
    # Exclude outliers (presumably due to measurement error)
    # rel water content should be between 0 to 100%
    filter(rel_water_content < 1.25) %>%
    filter(rel_water_content > -0.25) %>%
    select(-weight_pre, -weight_dry)
}

#' Calculate relative water content for sporophytes at the species level
#'
#' @param data Filmy fern relative water content data at the individual level
#' from desiccation tolerance test
#'
#' @return Dataframe
#' 
calculate_mean_water <- function (data) {
  
  data %>%
    group_by(salt, species, rec_time, dry_time) %>%
    summarize(
      rel_water = mean(rel_water_content),
      sd = sd(rel_water_content),
      .groups = "drop"
    )
  
}

#' Make a dataframe combining range size, vpd, and recovery
#' 
#' For preparing a dataframe for analysis with caper::pgls()
#'
#' @param combined_species_means Dataframe with desiccation tolerance data
#' @param env_range_data Dataframe with environmental (VPD) and range size data
#'
#' @return Datafrmae
#'
combine_env_env_range_recover <- function (combined_species_means, env_range_data) {

  combined_species_means %>%
  filter(!is.na(recovery_mean)) %>%
  select(species, generation, recovery = recovery_mean) %>%
  pivot_wider(names_from = "generation", values_from = "recovery", names_prefix = "recovery_") %>%
  left_join(
    select(env_range_data, 
           species, 
           sporo_vpd, 
           gameto_vpd,
           sporo_range_breadth,
           gameto_range_breadth
    ), 
    by = "species") %>%
  as.data.frame()
}

# Convert genus to just first letter
abbrev_sp <- function(data) {
  data %>%
    separate(species, c("genus", "epithet")) %>%
    mutate(genus = substr(genus, 1, 1)) %>%
    unite("species", genus, epithet)
}

# t-test ----

#' Run a t-test comparing response values between gametophytes and sporophytes.
#' 
#' Helper function for run_t_test()
#'
#' @param data Data with response variable (DT recovery, ETRmax or PAR95) for
#' filmy fern sporophytes and gametophytes in wide format. Three columns: `species`,
#' `gameto`, and `sporo`. `gameto` and `sporo` are each list-columns, where each
#' element of the list is a character vector of response variables.
#'
#' @return The results of running a two-sided t-test comparing response variables
#' between gametophytes and sporophytes.
#' 
run_t_test_gs <- function (data) {
  
  data %>%
    # Count number of samples per species per generation
    mutate(
      n_gameto = map_dbl(gametophyte, length),
      n_sporo = map_dbl(sporophyte, length)
    ) %>%
    # Only keep those with multiple samples
    filter(n_gameto > 1, n_sporo > 1) %>%
    # Run two-sided t-test, looping over species
    # see https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
    # for more info
    nest(data = c(gametophyte, sporophyte)) %>%
    mutate(
      t_test_res = map(data, ~t.test(x = .$sporophyte[[1]], y = .$gametophyte[[1]], paired = FALSE)),
      tidy_res = map(t_test_res, broom::tidy)
    ) %>%
    select(species, tidy_res) %>%
    unnest(cols = c(tidy_res)) %>%
    janitor::clean_names()
  
}

#' Run a t-test comparing response values between gametophytes and sporophytes
#' for DT recovery, PAR95, and ETRmax
#'
#' @param recovery_data Data with recovery from desiccation treatment with
#' one value per individual. Includes columns `species`, `salt`, `dry_time`, `rec_time`,
#' `recover`, and `generation`.
#' @param par_indiv Data with critical PAR (PAR where reach 95% of max ETR). One row
#' per individual. Includes columns `species`, `generation`, and `par_critical`
#' @param etr_indiv Data with maximum ETR. One row per individual. Includes columns
#' `species`, `generation`, and `etr_max`.
#'
#' @return Tibble. Columns include:
#' - `estimate`: estimate of difference between means (sporo value - gameto value)
#' - `estimate1`: estimate of mean value of sporophytes
#' - `estimate2`: estimate of mean value of gametophytes
#' - `statistic`: t-value
#' - `p_value`: p-value for the t-test
#' - `response`: response type ("recover", "par", or "etr")
#' 
run_t_test <- function (recovery_data, par_indiv, etr_indiv) {
  
  # Run t-test for DT recovery between sporophytes and gametophytes
  dt_t_test <-
    recovery_data %>%
    # Only compare recovery after 48hr between sporos and gametos dried using MgNO3 
    filter(salt == "MgNO3", rec_time == "48hr", dry_time == "2") %>%
    select(species, recovery, generation) %>%
    # Remove NA values
    filter(!is.na(recovery), !is.infinite(recovery)) %>%
    # Convert to wide format. Values in columns are now a list of numeric vectors
    # (recovery values)
    pivot_wider(names_from = "generation", values_from = "recovery", values_fn = list) %>% 
    run_t_test_gs() %>%
    mutate(response = "recovery")
  
  # Run t-test for PAR95 between sporophytes and gametophytes
  par_t_test <-
    par_indiv %>%
    select(species, generation, par_critical) %>%
    pivot_wider(names_from = "generation", values_from = "par_critical", values_fn = list) %>%
    run_t_test_gs() %>%
    mutate(response = "par")
  
  # Run t-test for ETRmax between sporophytes and gametophytes
  etr_t_test <-
    etr_indiv %>%
    select(species, generation, etr_max) %>%
    pivot_wider(names_from = "generation", values_from = "etr_max", values_fn = list) %>%
    run_t_test_gs() %>%
    mutate(response = "etr")
  
  # Combine the results
  bind_rows(
    dt_t_test,
    par_t_test,
    etr_t_test
  )
  
}

# Phylogentic signal ----

#' Analyze phylogenetic signal in a continuous trait of interest
#'
#' @param selected_trait Name of trait to analyze phylogenetic signal
#' @param traits Dataframe including all untransformed traits, with
#' 'species' as a column and other traits as other columns.
#' @param phy Phylogeny
#'
#' @return List of estimated Blomberg's K and Pagel's lambda and
#' their significance
#' 
analyze_cont_phylosig <- function (selected_trait, traits, phy) {
  
  # Trim data to non-missing trait values and
  # make sure species in same order in tree and traits
  traits_select <- traits %>% select(species, all_of(selected_trait)) %>%
    remove_missing(na.rm = TRUE)
  
  traits_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "traits") 
  phy_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "tree") 
  
  # Extract named vector of trait values for phylosig()
  trait_vec <- pull(traits_trim, selected_trait) %>%
    set_names(traits_trim$species)
  
  # Run phylosig() on selected trait
  # using Blomberg's K
  k_model <- phytools::phylosig(phy_trim, trait_vec, method = "K", test = TRUE)
  # and Pagel's lambda
  lambda_model <- phytools::phylosig(phy_trim, trait_vec, method = "lambda", test = TRUE)
  
  # get model results
  list(trait = selected_trait,
       kval = k_model$K,
       k_pval = k_model$P,
       lambda = lambda_model$lambda,
       lambda_pval = lambda_model$P)
  
}

#' Analyze phylogenetic signal by generation (sporophyte or gametophyte)
#'
#' @param combined_species_means Combined species means, including a 'species' column, columns for each trait,
#' and 'generation' column (either 'sporo' or 'gameto')
#' @param phy Phylogeny
#' @param traits_select Character vector of traits to analyze (must match column names in combined_species_means)
#'
#' @return Dataframe
#' 
analyze_phylosig_by_generation <- function(combined_species_means, phy, traits_select = c("recovery_mean", "etr_mean", "par_mean")) {
  
  # Subset to only gametophytes
  gameto_traits <- combined_species_means %>%
    filter(generation == "gametophyte")
  
  # Analyze phylogenetic signal
  gameto_phylosig <- map_df(traits_select, ~analyze_cont_phylosig(., gameto_traits, phy)) %>%
    mutate(generation = "gametophyte")
  
  # Subset to only sporophytes
  sporo_traits <- combined_species_means %>%
    filter(generation == "sporophyte")
  
  # Analyze phylogenetic signal
  sporo_phylosig <- map_df(traits_select, ~analyze_cont_phylosig(., sporo_traits, phy)) %>%
    mutate(generation = "sporophyte")
  
  # Combine results
  bind_rows(sporo_phylosig, gameto_phylosig)
  
}

# Modeling ----

#' Run (phylogenetic) generalized linear mixed models for filmy fern data
#' 
#' A phylogenetic effect is included for recovery from desiccation,
#' but not `etr` or `par`.
#'
#' @param combined_species_means Means calculated at the species level for response
#' variables ('recover_mean', 'etr_mean', and 'par_mean'), also including columns
#' for 'species' and 'generation'.
#' @param traits Trait data at the species level
#' @param phy Phylogenetic tree; list of class 'DNAbin'
#'
#' @return Dataframe; one row per model, with model testing for effect of
#' generation, growth habit, both (generation + growth habit), and their
#' interaction (generation * growth habit), or null
#' (effect of species only) on each response variable.
#' 
run_glmm <- function(combined_species_means, traits, phy) {
  
  # format trait data ---
  
  # - first combine DT and light response data 
  # and set habit, generation, range, and species as factors
  dt_light_data <-
    select(combined_species_means, species, generation, recovery = recovery_mean, etr = etr_mean, par = par_mean) %>%
    left_join(traits, by = "species") %>%
    mutate(
      habit = as.factor(habit),
      generation = as.factor(generation),
      species = as.factor(species))
  
  # then split each into a separate dataset with zero missing data
  # and drop unused levels
  dt_data <- dt_light_data %>%
    select(species, generation, habit, recovery) %>%
    remove_missing(na.rm = TRUE) %>%
    mutate_if(is.factor, droplevels)
  
  etr_data <- dt_light_data %>%
    select(species, generation, habit, etr) %>%
    remove_missing(na.rm = TRUE) %>%
    mutate_if(is.factor, droplevels)
  
  par_data <- dt_light_data %>%
    select(species, generation, habit, par) %>%
    remove_missing(na.rm = TRUE) %>%
    mutate_if(is.factor, droplevels)
  
  # finally, put these into a tibble for looping later over mcmcGLMM
  data_tibble <- tibble(
    response = c("etr", "par", "recovery"),
    data = list(etr_data, par_data, dt_data))
  
  # format phylogeny (only for DT) ---
  
  # have to remove node labels or inverseA won't work
  phy$node.label <- NULL
  
  # restrict phylogeny to only those species with DT data
  dt_phy <- ape::keep.tip(phy, as.character(dt_data$species))
  
  # set up inverse phylogenetic distance matrix
  dt_phy_inv <- MCMCglmm::inverseA(dt_phy, nodes="TIPS", scale = TRUE)
  
  # Run GLMMM ---
  
  # set number of iterations
  num <- 500000
  
  # define prior list. Need one G for each random effect, one R for each fixed effect
  my_prior <- list(G=list(G1=list(V=1,nu=0.02)),
                   R=list(R1=list(V=1,nu=0.02),R2=list(V=1,nu=0.02),R3=list(V=1,nu=0.02)))
  
  # run MCMCglmm for each dataset: non-phylogenetic models ---
  non_phylo_models <-
    list(
      fixed_effects = c("habit", "generation", "both", "intersect", "null"),
      response = c("etr", "par")
    ) %>%
    cross_df() %>%
    mutate(formula = case_when(
      fixed_effects == "null" ~ glue("{response}~1"),
      fixed_effects == "both" ~ glue("{response}~habit+generation"),
      fixed_effects == "intersect" ~ glue("{response}~habit*generation"),
      fixed_effects == "habit" ~ glue("{response}~habit"),
      fixed_effects == "generation" ~ glue("{response}~generation")
    )) %>%
    mutate(formula = map(formula, as.formula)) %>%
    left_join(data_tibble, by = "response") %>%
    mutate(
      model = map2(
        formula, 
        data,
        ~MCMCglmm::MCMCglmm(
          fixed = .x, 
          random = ~species, 
          family = "gaussian", 
          prior = my_prior, 
          data = as.data.frame(.y),
          nitt = num, burnin = 1000, thin = 500,
          verbose = FALSE)
      )
    )
  
  # run MCMCglmm for each dataset: phylogenetic models ---
  phylo_models <-
    list(
      fixed_effects = c("habit", "generation", "both", "intersect", "null"),
      response = c("recovery")
    ) %>%
    cross_df() %>%
    mutate(formula = case_when(
      fixed_effects == "null" ~ glue("{response}~1"),
      fixed_effects == "both" ~ glue("{response}~habit+generation"),
      fixed_effects == "intersect" ~ glue("{response}~habit*generation"),
      fixed_effects == "habit" ~ glue("{response}~habit"),
      fixed_effects == "generation" ~ glue("{response}~generation")
    )) %>%
    mutate(formula = map(formula, as.formula)) %>%
    left_join(data_tibble, by = "response") %>%
    mutate(
      model = map2(
        formula, 
        data,
        ~MCMCglmm::MCMCglmm(
          fixed = .x, 
          random = ~species, 
          family = "gaussian", 
          ginverse = list(species = dt_phy_inv$Ainv), 
          prior = my_prior, 
          data = as.data.frame(.y),
          nitt = num, burnin = 1000, thin = 500,
          verbose = FALSE)
      )
    )
  
  bind_rows(non_phylo_models, phylo_models)
  
}

# define function to extract p value table from summary as dataframe, append model name
summary_to_csv <- function (model) {
  results.df <- as.data.frame(summary(model)$solutions)
  results.df$model <- deparse(substitute(model))
  results.df$effect <- rownames(results.df)
  rownames(results.df) <- NULL
  colnames(results.df) <- c("Parameter estimate", "Lower 95% CI", "Upper 95% CI", "Effective sample size", "P-value", "Model", "Effect")
  return(results.df)
}

#' Run phylogenetic generalized least squares (PGLS)
#' for range size and VPD vs. desiccation tolerance
#' in filmy ferns from Moorea
#'
#' @param env_range_recover_data Dataframe with range size, environmental (VPD), 
#' and desiccation tolerance data
#' @param phy Phylogeny
#'
#' @return list of model objects:
#' - sporo_vpd_model:
#' - gameto_vpd_model:
#' - sporo_range_model:
#' - gameto_range_model:
#' 
run_pgls <- function (env_range_recover_data, phy) {
  
  # Create a comparative.data object for caper combining the range size, vpd, 
  # and recovery data with the phylogeny
  env_range_recover_data_comp <- caper::comparative.data(
    phy, env_range_recover_data, species, na.omit = FALSE)
  
  # Run PGLS for each combination of generation x vpd or range breadth vs. recovery
  res <- list(
    sporo_vpd_model = caper::pgls(recovery_sporophyte ~ sporo_vpd, env_range_recover_data_comp),
    gameto_vpd_model = caper::pgls(recovery_gametophyte ~ gameto_vpd, env_range_recover_data_comp),
    sporo_range_model = caper::pgls(recovery_sporophyte ~ sporo_range_breadth, env_range_recover_data_comp),
    gameto_range_model = caper::pgls(recovery_gametophyte ~ gameto_range_breadth, env_range_recover_data_comp)
    )
  
}

#' Extract fitted values from a PGLS model
#'
#' @param model The PGLS model (list of class "pgls")
#'
#' @return Tibble, with two columns: species and the predicted value
#' 
predict_pgls <- function(model) {
  model$fitted %>%
    as.data.frame() %>%
    rownames_to_column("species") %>%
    rename(predicted = V1) %>%
    as_tibble()
}

#' Summarize statistics of a PGLS model
#'
#' @param model The PGLS model (list of class "pgls")
#'
#' @return Tibble, with various columns for model statistics (r_squared, p-value, etc)
#' 
tidy_pgls <- function (model) {
  model_sum <- summary(model)
  tibble(
    sigma = model_sum$sigma,
    df = model_sum$df %>% paste(collapse = ","),
    r_squared = model_sum$r.squared,
    f_value = model_sum[["fstatistic"]][["value"]],
    adj_r_squared = model_sum$adj.r.squared[,1],
    # Need to run anova() to get p-value
    p_value = anova(model) %>% filter(!is.na(`F value`)) %>% pull(`Pr(>F)`)
  )
}

# Etc ----

#' Match trait data and tree
#' 
#' Order of species in traits will be rearranged to match the
#' phylogeny.
#'
#' @param traits Dataframe of traits, with 'species' column and
#' additional columns, one for each trait
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the traits, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
match_traits_and_tree <- function (traits, phy, return = c("traits", "tree")) {
  
  assert_that("species" %in% colnames(traits))
  
  # Keep only species in phylogeny
  traits <- traits %>%
    filter(species %in% phy$tip.label) 
  
  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, traits$species))
  
  # Get traits in same order as tips
  traits <- left_join(
    tibble(species = phy$tip.label),
    traits,
    by = "species"
  )
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(traits$species, phy$tip.label)))
  
  # Return traits or tree
  assert_that(return %in% c("tree", "traits"))
  
  if(return == "tree") { 
    return (phy) 
  } else {
    return (traits)
  }
  
}

#' Match community data and tree
#' 
#' Order of species in comm will be rearranged to match the
#' phylogeny.
#'
#' @param comm Community data frame, with one column for sites and
#' the rest for species.
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the community, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
match_comm_and_tree <- function (comm, phy, return = c("comm", "tree")) {
  
  assert_that("species" %in% colnames(comm))
  
  # Keep only species in phylogeny
  comm <- comm %>%
    filter(species %in% phy$tip.label) 
  
  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, comm$species))
  
  # Get comm in same order as tips
  comm <- left_join(
    tibble(species = phy$tip.label),
    comm
  )
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, phy$tip.label)))
  
  # Return comm or tree
  assert_that(return %in% c("tree", "comm"))
  
  if(return == "tree") { 
    return (phy) 
  } else {
    return (comm)
  }
  
}

# Plots ----

# Format legend (include title, but get rid of boxes)
legend_theme <- theme (legend.background = element_rect(colour = "transparent", fill = "transparent", size = 0.5),
                       legend.key = element_rect(colour = "white", fill = "white", size = 0.5),
                       legend.text = element_text(size = 12))

make_sporo_dt_plot <- function  (recovery_species_means, traits) {
  # Define color-blind palette
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  plot_data <-
    recovery_species_means %>%
    # Use only sporophytes
    filter(generation == "sporo") %>%
    # Don't plot control
    filter(salt != "Control") %>%
    # Add growth habit data
    left_join(traits, by = "species") %>%
    assert(not_na, habit, range) %>%
    # Abbreviate species names to G. species
    mutate(species = str_replace_all(species, "[a-z]+_", ". ")) %>%
    # Reorder levels for species: alphabetical within growth habit
    mutate(species = fct_reorder2(species, habit, species, .desc = FALSE)) %>%
    rename(recover = recover_mean, sd = recover_sd, n = recover_n)
  
  # Self-defined formatting function for percent
  percent_formatter <- function(x) {
    lab <- x*100
  }

  # Set colors for desiccation intensity
  col <- c(cbPalette[7], cbPalette[5], cbPalette[6])
  
  # Set offset of datapoints (dodge)
  pd <- position_dodge(.3)
  
  # Make plot
  ggplot(data = plot_data, aes(x = rec_time, y = recover, group = interaction(dry_time, salt), shape = salt, fill = salt)) +
    geom_line(position = pd, aes(linetype = dry_time)) +
    geom_errorbar(
      aes(ymin = recover - sd, ymax = recover + sd), 
      colour = "grey50", width = 0, 
      position = pd) +
    geom_point(position = pd, size = 1.5) +
    scale_linetype_manual(name = "Desiccation Time",
                          values = c("solid", "dashed"),
                          breaks = c("15d", "2d"),
                          labels = c("15 days", "2 days")) +
    scale_shape_manual(name = "Desiccation Intensity",
                       values = c(23,24,21),
                       breaks = c("NaCl", "MgNO3",  "LiCl"),
                       labels = c("-38 MPa", "-86 MPa", "-282 Mpa")) + 
    scale_fill_manual(name = "Desiccation Intensity",
                      values = col, 
                      breaks = c("NaCl", "MgNO3",  "LiCl"),
                      labels = c("-38 MPa", "-86 MPa", "-282 Mpa")) + 
    scale_y_continuous(label = percent_formatter, limits = c(0,1.1), expand = c(0,0), breaks=c(0,.2,.4,.6,.8,1.0)) +
    scale_x_discrete(limits = c("30min","24hr","48hr"),
                     labels = c("0.5", "24", "48")) +
    xlab("Recovery Time (hr)") +
    ylab("Recovery (%)") +
    theme_bw() +
    legend_theme +
    theme(legend.position="bottom",
          panel.grid.minor=element_blank(), 
          panel.grid.major=element_blank()) + 
    facet_wrap( ~ species, ncol=3)
}

#' Determine which species have widespread gametophytes
#'
#' @param raw_community_matrix Community matrix of Nitta et al 2017
#' @param filmy_species Character vector of filmy fern species names
#' @param phy Phylogeny including filmy ferns of Moorea and Tahiti
#' @param moorea_sites Site data for plots on Moorea including elevation
#' 
#' @return Dataframe
#' 
analyze_dist_pattern <- function(community_matrix_raw, filmy_species, phy, moorea_sites) {
  
  # Convert raw community data to long format by sporophyte and gametophyte
  # min and max elevational range.
  range_long <-
    community_matrix_raw %>%
    rename(site = X1) %>% 
    gather(species, abundance, -site) %>%
    filter(species %in% filmy_species) %>%
    mutate(generation = case_when(
      str_detect(site, "_S") ~ "sporophyte",
      str_detect(site, "_G") ~ "gametophyte"
    )) %>%
    mutate(site = str_remove_all(site, "_S") %>% str_remove_all("_G")) %>%
    # Keep only sites on Moorea
    filter(site %in% moorea_sites$site) %>%
    left_join(moorea_sites, by = "site") %>%
    filter(abundance > 0) %>%
    assert(not_na, el) %>%
    group_by(species, generation) %>%
    summarize(
      min_range = min(el),
      max_range = max(el),
      .groups = "drop"
    )
  
  # Compare ranges to see which species have gametophytes beyond sporophytes
  full_join(
    range_long %>%
      gather(variable, value, -species, -generation) %>%
      filter(generation == "sporophyte") %>% 
      select(-generation) %>%
      spread(variable, value) %>%
      rename_at(vars(max_range, min_range), ~paste("sporo", ., sep = "_")),
    range_long %>%
      gather(variable, value, -species, -generation) %>%
      filter(generation == "gametophyte") %>% 
      select(-generation) %>%
      spread(variable, value) %>%
      rename_at(vars(max_range, min_range), ~paste("gameto", ., sep = "_")),
    by = "species"
  ) %>%
    assert(is_uniq, species) %>%
    mutate(
      range = case_when(
        gameto_max_range > sporo_max_range + 200 ~ "widespread",
        gameto_min_range < sporo_min_range - 200 ~ "widespread",
        is.na(sporo_max_range) & is.na(sporo_min_range) ~ "widespread",
        TRUE ~ "not_widespread"
      )
    ) %>%
    mutate(gameto_range_breadth = gameto_max_range - gameto_min_range) %>%
    mutate(sporo_range_breadth = sporo_max_range - sporo_min_range)
  
}

# Manuscript ----

# Helper function to check on p-value significance and extract from t_test_results
#
#' @param species_select Name of species to select from t-test results
#' @param response_select Name of response to select from t-test results
#' @param data t-test results
#' @param p_check Boolean: should the p-value be required to be less than 0.05?
#'
#' @return A p-value for the selected species and response
ttest_pval <- function(species_select, response_select, data = t_test_results, p_check = TRUE) {
  data %>% 
    # Filter to appropriate species and response
    filter(species == species_select, response == response_select) %>%
    # Verify the p-value is significant (or override with `p_check = FALSE`)
    assert(function (x) magrittr::is_less_than(x, 0.05) | !p_check, p_value) %>%
    # Format the p-value
    mutate(
      p_value = case_when(
        p_value < 0.001 ~ "<0.001",
        TRUE ~ jntools::round_t(p_value, digits = 3))
    ) %>% 
    pull(p_value)
}

# Dummy function to track arbitary output from rmarkdown::render()
render_tracked <- function (tracked_output, dep1 = NULL, dep2 = NULL, ...) {
  rmarkdown::render(...)
}

#' Convert latex file to docx
#' 
#' Requires pandoc to be installed and on command line
#'
#' @param latex Path to input latex file.
#' @param docx Path to output docx file.
#' @param template Optional; path to template docx file.
#' @param lua_filter Optional; path to lua filter file.
#' @param wd Working directory to run conversion. Should be same as
#' directory containing any files needed to render latex to pdf.
#' @param ... Other arguments; not used by this function, but for tracking dependencies.
#'
#' @return List including STDOUT of pandoc; externally, the
#' docx file will be rendered in `wd`.
#' 
latex2docx <- function (latex, docx, template = NULL, lua_filter = NULL, wd = getwd(), ...) {
  
  assertthat::assert_that(assertthat::is.readable(latex))
  
  assertthat::assert_that(assertthat::is.dir(fs::path_dir(docx)))
  
  latex <- fs::path_abs(latex)
  
  docx <- fs::path_abs(docx)
  
  template <- if (!is.null(template)) {
    glue::glue("--reference-doc={fs::path_abs(template)}")
  } else {
    NULL
  }
  
  lua_filter <- if (!is.null(lua_filter)) {
    glue::glue("--lua-filter={fs::path_abs(lua_filter)}")
  } else {
    NULL
  }
  
  processx::run(
    command = "pandoc",
    args = c("-s", latex, template, "-o", docx, lua_filter),
    wd = wd
  )
  
}

