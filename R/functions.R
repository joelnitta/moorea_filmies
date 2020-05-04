# Define functions used in the analysis

clean_recovery <- function (recovery_data_raw) {
  
  # Convert any recovery >100% to 1
  recovery_data_raw %>%
    mutate(recover = case_when(
      recover > 1 ~ 1,
      TRUE ~ recover
    ))
  
}

summarize_recovery <- function (recovery_data) {
  
  # Summarize recovery by salt treatment, dry time, recovery time, generation, and species
  recovery_data %>%
    group_by(salt, drytime, rectime, generation, genus_species) %>%
    summarize(
      mean_recover = mean(recover, na.rm = TRUE),
      sd_recover = sd(recover, na.rm = TRUE),
      n = n()
    )
  
}


# filmy_process_physio
# Calculates mean values per species from raw physiological data 
# (percent desiccation tolerance and light curve data)

# input: raw DT and light response data
# output: list of subsetted raw DT data, mean DT by species, 
# ETR by individual/species, and PAR by individual/species

# uses summarySE, which gets loaded by script that calls process_physio

process_physio <- function (DT_data, light_data) {
  
  ######################################################
  ### Desiccation Tolerance (48 recovery from MGNO3) ###
  ######################################################
  
  # use temporary "data" object to store raw DT data
  data <- DT_data
  
  ### miscellaneous data fixes
  
  # for sake of analysis, convert Calllistopteris apiifolia NaCl treatment to MGNO3
  Capi <- data[data$genus_species == "Callistopteris_apiifolia" & data$generation == "sporo" & data$salt == "NaCl",]
  Capi$salt <- "MgNO3"
  data <- rbind(data, Capi)
  rm(Capi)
  
  # for sake of analysis, convert Abrodictyum dentatum NaCl treatment to MGNO3
  Aden <- data[data$genus_species == "Abrodictyum_dentatum" & data$generation == "sporo" & data$salt == "NaCl",]
  Aden$salt <- "MgNO3"
  data <- rbind(data, Aden)
  rm(Aden)
  
  # convert any recovery >100% to 1
  data$recover[data$recover>1] <- 1
  
  # convert from digit to percent
  data$recover <- data$recover * 100
  
  ### subset data: compare filmy gameto and sporo, MgNO3 2d treatment
  data <- data[data$salt == "MgNO3",]
  data <- data[data$drytime == "2d",]
  data <- data[data$rectime == "48hr",]
  
  # drop unused levels
  data <- droplevels(data)
  
  ### use summarySE function to get mean and SE by species and generation
  mean.DT <- summarySE(data, measurevar="recover", groupvars=c("genus_species", "generation"))
  
  # keep subsetted raw DT data
  raw.DT <- data
  
  # going to use "data" again for light dataset, so get rid of DT data
  rm(data)
  
  ###############
  ### Max ETR ###
  ###############
  
  # use temporary "data" object to store raw light data
  data <- light_data
  
  # subset data (drop outlier points / samples)
  data <- data[data$exclude_sample == 0,]
  data <- data[data$exclude_points == 0,]
  
  # drop unused levels
  data <- droplevels(data)
  
  ### set up loop ###
  
  # use split function to split raw data frame into list of smaller dataframes by individual
  df_list <- split(data, as.factor(data$individual))
  
  ### calculate ETR max ###
  # loop to get max ETR for each dataframe in the list
  genus_species.loop = NULL
  condition.loop = NULL
  generation.loop = NULL
  out.ETR = NULL
  individual.loop = NULL
  
  for (i in 1:length(df_list)) { 
    # extract maximum ETR
    out.ETR[i] <- max(df_list[[i]]$ETR)
    # include individual
    individual.loop[i] <- as.character(unique(df_list[[i]]$individual))
    # include genus and species
    genus_species.loop[i] <- as.character(unique(df_list[[i]]$genus_species))
    # include gameto/sporo label
    generation.loop[i] <- as.character(unique(df_list[[i]]$generation))
  }
  
  # summary of ETR by individual
  summary.ETR <- data.frame("individual" = individual.loop, "ETRmax" = out.ETR, "genus_species" = genus_species.loop, "generation" = generation.loop)
  
  # summary of mean values by species.
  mean.ETR <- summarySE(summary.ETR, measurevar="ETRmax", groupvars=c("genus_species", "generation"))
  
  ### PAR95 ###
  
  # use same df_list that we created for ETR analysis above
  
  # set up loop
  # for dummy variables that are lists, need to create vector of lists, with length equal to total number of repeats
  data.sub = vector('list',length(df_list))
  out = vector('list',length(df_list))
  out.sum = vector('list',length(df_list))
  k.est = NULL
  PAR.crit = NULL
  filename = NULL
  fitline = NULL
  individual.loop = NULL
  genus_species.loop = NULL
  out.ETR = NULL
  generation.loop = NULL
  
  # loop to get critical PAR value for each dataframe in the list
  
  # Specify the function as ETR=A(1-exp(-k*PAR))
  # Choose A to be equal to the maximum ETR value
  # Then do non-linear least squares on k
  
  for (i in 1:length(df_list)) { 
    # nonlinear least-squares estimates of parameters of the subsetted data
    out[[i]] <- nls(ETR~max(ETR)*(1-exp(-k*PAR)),start=list(k=0.04),data=df_list[[i]],trace=F,control = list(maxiter = 500))
    # summarize results
    out.sum[[i]] <- summary(out[[i]])
    # grab estimated value of k
    k.est[i] <- out.sum[[i]]$parameters[, "Estimate"]
    # calculate critical PAR (PAR where reach 95% of max ETR)
    PAR.crit[i] <- -log(0.05)/k.est[i]
    # include individual
    individual.loop[i] <- as.character(unique(df_list[[i]]$individual))
    # include genus and species
    genus_species.loop[i] <- as.character(unique(df_list[[i]]$genus_species))
    # include gameto/sporo label
    generation.loop[i] <- as.character(unique(df_list[[i]]$generation))
  }
  
  # summary by individuals
  summary.PAR <- data.frame("individual" = individual.loop, "PAR95" = PAR.crit, "K" = k.est, "genus_species" = genus_species.loop, "generation" = generation.loop)
  
  # summary of mean values by species.
  mean.PAR <- summarySE(summary.PAR, measurevar="PAR95", groupvars=c("genus_species", "generation"))
  
  results <- list (raw.DT, mean.DT, summary.ETR, mean.ETR, summary.PAR, mean.PAR)
  names(results) <- c("raw.DT", "mean.DT", "summary.ETR", "mean.ETR", "summary.PAR", "mean.PAR" )
  return(results)
}

## Function "Summary SE" summarizes data.
## Downloaded from: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm),
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean"=measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
