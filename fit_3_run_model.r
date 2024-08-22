# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              
#                 Forest Regrowth Model Fitting and Comparison                 
#                                                                              
#                            Ana Avila - August 2024                           
#                                                                              
#     This script runs and compares the Chapman-Richards growth curve fit      
#     through optim, as well as linear, GAM, and random forest models for      
#     forest regrowth analysis. The results are combined into a single         
#     dataframe and saved as 'region_results_all.csv' in the            
#     './data/' directory.    
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)

# Source external R scripts for data import and function definitions
source("fit_1_import_data.r")
source("fit_2_functions.r")

set.seed(1)
ncores <- 30
registerDoParallel(cores = ncores)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------------- Switches ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

test_switch <- FALSE # Set to TRUE for test_rows or FALSE to run all rows
# include test_rows rows to test
test_rows <- 1
split_biome <- FALSE # FALSE to execute for the whole country, TRUE to split dataframe between biomes

fit_logistic <- TRUE

run_initial_fit_by_order <- TRUE

optim_switch <- TRUE
lm_switch <- FALSE
rf_switch <- FALSE
all_switch <- FALSE

# region <- "countrywide"
region <- "amaz"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# number of rows to be included in analysis
n_samples <- 10000

# name for export
name <- region
if (split_biome) {
  name <- paste0(name, "_split_biome")
}


# Define climatic parameters that change yearly
climatic_pars <- c("prec", "si")

# Define parameters that do not correspond to data, used for functional form
if (fit_logistic) {
  name <- paste0(name, "_logistic")
  conditions <- list()
  non_data_pars <- c("k0", "B0")
} else {
  # Define conditions for parameter constraints
  conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0')
  non_data_pars <- c("k0", "B0_exp", "B0", "theta")
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define land-use history intervals to import four dataframes
intervals <- list("5yr", "10yr", "15yr", "all")
datafiles <- paste0("./data/", region, "_", intervals, ".csv")
dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparing dataframes whether or not there is a split between biomes

if (split_biome) {
  # 1 = Amazonia
  # 2 = Caatinga
  # 3 = Cerrado
  # 4 = Mata Atlantica
  # 5 = Pampa
  # 6 = Pantanal
  biomes <- c("amaz", "atla")

  dataframes <- lapply(dataframes, function(df) {
    df[df$biome %in% c(1, 4), ]
  })

  # Split each dataframe by biome and store in a nested list
  dataframes <- lapply(dataframes, function(df) {
    split(df, df$biome)
  })

  # Assuming value_sample is defined as the number of rows to sample
  dataframes <- lapply(dataframes, function(list_of_dfs) {
    lapply(list_of_dfs, function(df) {
        df_sampled <- df[sample(nrow(df), n_samples, replace = FALSE), ]
      return(df_sampled)
    })
  })

} else {
  dataframes <- lapply(dataframes, function(df) {
    df_sampled <- df[sample(nrow(df), n_samples, replace = FALSE), ]
    return(df_sampled)
  })
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Identify common columns across all dataframes
# (avoids errors for parameters like last_LU, that not all dataframes may contain depending on how many rows were sampled)
# and filter out non-predictors

if (split_biome) {

  colnames_lists <- lapply(dataframes, function(df_list) {
    lapply(df_list, function(df) colnames(df))
  })

  colnames_lists <- lapply(colnames_lists, unlist)

  # Step 2: Find the intersection of all column names
  colnames_intersect <- Reduce(intersect, colnames_lists)

} else {
  colnames_intersect <- Reduce(intersect, map(dataframes, colnames))
}


colnames_filtered <- colnames_intersect[!grepl(
  "age|agbd|latitude|longitude|prec|si|mature_biomass|distance|biome|cwd",
  colnames_intersect
)]

# Define different sets of data parameters for modeling
data_pars <- list(
  c("cwd"),
  c("cwd", "mean_prec", "mean_si"),
  c("cwd", climatic_pars),
  colnames_filtered[!grepl("ecoreg|soil", colnames_filtered)], # land use only
  colnames_filtered,
  c(colnames_filtered, "cwd", "mean_prec", "mean_si"),
  c(colnames_filtered, "cwd", climatic_pars)
)

data_pars_names <- c(
  "cwd", "meanclim", "yearlyclim", "LU", "LU_soil_ecor",
  "LU_soil_ecor_meanclim", "LU_soil_ecor_yearlyclim"
)

# Define basic parameter sets for modeling
basic_pars <- list(
  c("age", "B0"),
  c("age", "k0", "B0"),
  c("k0", "B0"), # age is implicit
  c("k0", "B0_exp") # age is implicit
)

if (fit_logistic) {
  basic_pars <- basic_pars[1:3]
}

basic_pars_names <- as.list(sapply(basic_pars, function(par_set) paste(par_set, collapse = "_")))

# for linear model and random forest
# - remove climatic_pars (since they are only used in optim for yearly change)
# - add mature_biomass (nearest neighbor) to each parameter set
# - add new parameter sets
data_pars_lm <- c(
  lapply(
    Filter(function(x) !any(climatic_pars %in% x), data_pars), # remove climatic_pars
    function(x) c(x, "mature_biomass")#, "age")
  ),
  list(
    c("age"), c("mature_biomass"), c("age", "mature_biomass")
  )
)

filtered_data_pars_names <- data_pars_names[!sapply(data_pars, function(x) any(climatic_pars %in% x))]
data_pars_lm_names <- c(filtered_data_pars_names, "age", "mat_biomass", "age_mat_biomass")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create grids of different combinations of model inputs

# Optim / growth curve
iterations_optim <- expand.grid(
  interval = seq_along(intervals),
  data_par = seq_along(data_pars),
  basic_par = seq_along(basic_pars)
)

basic_pars_with_age <- which(sapply(basic_pars, function(x) "age" %in% x))
data_pars_with_climatic <- which(sapply(data_pars, function(x) any(climatic_pars %in% x)))
# Remove rows where both conditions are met
iterations_optim <- iterations_optim %>% filter(!(basic_par %in% basic_pars_with_age & data_par %in% data_pars_with_climatic))

# Linear Model
iterations_lm <- expand.grid(
  interval = seq_along(intervals),
  data_par = seq_along(data_pars_lm)
)

if (split_biome) {
  # Create a grid for the biome dimension
  biome_grid <- expand.grid(
    biome = seq_along(biomes)
  )

  # Merge base_grid with biome_grid to include the biome dimension
  iterations_optim <- merge(iterations_optim, biome_grid, by = NULL)
  iterations_lm <- merge(iterations_lm, biome_grid, by = NULL)
}

# Random Forest
# Filter out single-predictor cases
rows_to_remove <- sapply(iterations_lm$data_par, function(i) {
  length(data_pars_lm[[i]]) == 1
})
iterations_rf <- iterations_lm[!rows_to_remove, ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------------------- Find ideal parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if (run_initial_fit_by_order){
  # keep combination of parameters for initial fit
  initial_pars <- find_combination_pars(iterations_optim, "nls")
} else {
  initial_pars <- readRDS(paste0("./data/", name, "_ideal_par_combination.rds"))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Run Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# test_rows <- which(iterations_optim$interval == 1 & iterations_optim$data_par == 5 &
#   iterations_optim$basic_par == 1)

start_time <- Sys.time()
print(start_time)

# Run linear model
if (lm_switch || all_switch) {
  results_lm <- run_foreach(iterations_lm, "lm", run_lm)
}

# Run optimization with growth curve
if (optim_switch || all_switch) {
  results_optim <- run_foreach(iterations_optim, "optim", run_optim, conditions)
}

# Run random forest
if (rf_switch || all_switch) {
  results_rf <- run_foreach(iterations_rf, "rf", run_rf)
}

# Combine all results if all_switch is TRUE
if (all_switch) {
  
  results_all <- bind_rows(results_optim, results_lm, results_rf) %>%
    arrange(data_pars, basic_pars, data_name) # sort rows by parameter

  
  write.csv(results_all, paste0("./data/", name, "_results_all.csv"), row.names = FALSE)

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

print(paste(
  region, "ALL DONE! Time for the whole thing using ", ncores, "cores: ",
  as.numeric(difftime(Sys.time(), start_time, units = "hours")), "hours"
))
