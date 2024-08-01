################################################################################
#                                                                              #
#                 Forest Regrowth Model Fitting and Comparison                 #
#                                                                              #
#                              Ana Avila - July 2024                           #
#                                                                              #
#     This script runs and compares the Chapman-Richards growth curve fit      #
#     through optim, as well as linear, GAM, and random forest models for      #
#     forest regrowth analysis. The results are combined into a single         #
#     dataframe and saved as 'fit_results_test_train_rf.csv' in the            #
#     './data/' directory.                                                     #
#                                                                              #
################################################################################


library(ggplot2)
library(foreach)
library(doParallel)
library(mgcv)
library(randomForest)

# Source external R scripts for data import and function definitions
source("fit_1_import_data.r")
source("fit_2_functions.r")

set.seed(1)
registerDoParallel(cores = 25)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----------------- Switches ---------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

test_switch <- FALSE # Set to TRUE for single row or FALSE for all rows
split_biome <- FALSE # FALSE to execute for the whole country, TRUE to split dataframe between biomes

optim_switch <- FALSE
lm_switch <- FALSE
rf_switch <- FALSE
all_switch <- FALSE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----------------- Global Variables ---------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define climatic parameters that change yearly
climatic_pars <- c("prec", "si")

# Define parameters that do not correspond to data, used for functional form
non_data_pars <- c("k0", "B0_exp", "B0", "theta")

# Define categorical parameters
pars_categ <- c("indig", "protec", names(data)[str_detect(names(data), "LU|ecoreg|soil")])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------------ Import Data -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define land-use history intervals to import four dataframes
intervals <- list("5y", "10y", "15y", "all")
datafiles <- paste0("./data/", intervals, "_LULC_countrywide.csv") # "_LULC_amaz.csv")
dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)
nrow(dataframes[[1]])

lapply(dataframes, nrow)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to create training and testing samples

sample_data <- function(df, size_train, size_test) {
  train_dataframes <- sample_n(df, size = size_train)
  remaining <- df[!(row.names(df) %in% row.names(train_dataframes)), ]
  test_dataframes <- sample_n(remaining, size = size_test)
  list(train_dataframes, test_dataframes)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparing dataframes whether or not there is a split between biomes

if (split_biome) {
  # 1 = Amazonia
  # 2 = Caatinga
  # 3 = Cerrado
  # 4 = Mata Atlantica
  # 5 = Pampa
  # 6 = Pantanal
  biomes <- c("amaz", "caat", "cerr", "mata", "pamp", "pant")

  # Split each dataframe by biome and store in a nested list
  split_dataframes <- lapply(dataframes, function(df) {
    split(df, df$biome)
  })

  # Initialize lists to store training and testing dataframes for each biome
  train_dataframes <- list()
  test_dataframes <- list()

  # Iterate over each biome and interval
  for (biome in biomes) {
    train_dataframes[[biome]] <- list()
    test_dataframes[[biome]] <- list()
    for (interval in seq_along(intervals)) {
      df_biome <- split_dataframes[[interval]][[biome]]
      samples <- sample_data(df_biome, size_train = 60000, size_test = 40000)
      train_dataframes[[biome]][[intervals[interval]]] <- samples[[1]]
      test_dataframes[[biome]][[intervals[interval]]] <- samples[[2]]
    }
  }
}

samples <- lapply(dataframes, sample_data, size_train = 60000, size_test = 40000)
train_dataframes <- lapply(samples, `[[`, 1)
test_dataframes <- lapply(samples, `[[`, 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------- Define Parameters ----------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define conditions for parameter constraints
conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0')

# Identify common columns across all dataframes
# (avoids errors for parameters like last_LU, that not all dataframes may contain depending on how many rows were sampled)
# and filter out non-predictors
colnames_intersect <- Reduce(intersect, map(dataframes, colnames))
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

basic_pars_names <- as.list(sapply(basic_pars, function(par_set) paste(par_set, collapse = "_")))

# for linear model and random forest
# - remove climatic_pars (since they are only used in optim for yearly change)
# - add mature_biomass (nearest neighbor) to each parameter set
# - add new parameter sets
data_pars_lm <- c(
  lapply(
    Filter(function(x) !any(climatic_pars %in% x), data_pars), # remove climatic_pars
    function(x) c(x, "mature_biomass")
  ),
  list(
    c("age"), c("mature_biomass"), c("age", "mature_biomass")
  )
)

filtered_data_pars_names <- data_pars_names[!sapply(data_pars, function(x) any(climatic_pars %in% x))]
data_pars_lm_names <- c(filtered_data_pars_names, "age", "mat_biomass", "age_mat_biomass")
print(data_pars_lm_names)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create grids of different combinations of model inputs

# Optim / growth curve
iterations_optim <- expand.grid(
  dataframe = seq_along(dataframes),
  data_par = seq_along(data_pars),
  basic_par = seq_along(basic_pars)
)

# Linear Model
iterations_lm <- expand.grid(
  dataframe = seq_along(dataframes),
  data_par = seq_along(data_pars_lm)
)

# Random Forest
# Filter out single-predictor cases
rows_to_remove <- sapply(iterations_lm$data_par, function(i) {
  length(data_pars_lm[[i]]) == 1
})
iterations_rf <- iterations_lm[!rows_to_remove, ]

print(iterations_rf)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------------- Run Model --------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

start_time <- Sys.time()
print(start_time)

prefix <- 'countrywide_'

# ~~~~~~~~~~~~~~~~ OPTIM ~~~~~~~~~~~~~~~~~~~~~#
# Run optimization with growth curve

if (optim_switch || all_switch) {

  results_optim <- foreach(iter = 1:nrow(iterations_optim), 
    .combine = "bind_rows", .packages = c("dplyr")) %dopar%
    {

    print(iter)
    i <- iterations_optim$dataframe[iter]
    j <- iterations_optim$data_par[iter]
    k <- iterations_optim$basic_par[iter]

    train_data <- train_dataframes[[i]]
    test_data <- test_dataframes[[i]]
    pars_chosen <- data_pars[[j]]
    pars_basic <- basic_pars[[k]]

    optim_output <- run_optim("nls", pars_basic, pars_chosen, train_data, test_data, conditions)

    # Organizes result into a new row for the final dataframe
    optim_row <- process_row(
      optim_output, "optim", intervals[[i]],
      data_pars_names[[j]],
      basic_pars_names = basic_pars_names[[k]]
    )

    print(optim_row)
    print(paste("Time so far: ", as.numeric(difftime(Sys.time(), start_time, units = "mins")), " minutes"))
    optim_row
  }

  df_optim <- write_results_to_csv(results_optim, prefix)
}




# ~~~~~~~~~~~~~~~~ LINEAR MODEL ~~~~~~~~~~~~~~~~~~#

if (lm_switch || all_switch) {

  results_lm <- foreach(iter = 1:nrow(iterations_lm), 
    .combine = "bind_rows", .packages = c("dplyr")
  ) %dopar% {

    print(iter)
    i <- iterations_lm$dataframe[iter]
    j <- iterations_lm$data_par[iter]

    train_data <- train_dataframes[[i]]
    test_data <- test_dataframes[[i]]
    pars_chosen <- data_pars_lm[[j]]

    lm_output <- run_lm(train_data, test_data, pars_chosen)

    # Organizes result into a new row for the final dataframe
    lm_row <- process_row(
      lm_output$pars, intervals[[i]],
      "lm", data_pars_lm_names[[j]], lm_output
    )

    print(lm_row)
    print(paste("Time so far: ", as.numeric(difftime(Sys.time(), start_time, units = "mins")), " minutes"))
    lm_row
  }

  df_lm <- write_results_to_csv(results_lm, prefix)
}

# ~~~~~~~~~~~~~~~~ RANDOM FOREST ~~~~~~~~~~~~~~~~~~#

if (rf_switch || all_switch) {

  results_rf <- foreach(iter = 1:nrow(iterations_rf), 
    .combine = "bind_rows", .packages = c("dplyr", "randomForest")
  ) %dopar% {

    print(iter)
    i <- iterations_rf$dataframe[iter]
    j <- iterations_rf$data_par[iter]

    train_data <- train_dataframes[[i]]
    test_data <- test_dataframes[[i]]
    pars_chosen <- data_pars_lm[[j]]

    rf_output <- run_rf(train_data, test_data, pars_chosen)

    # Organizes result into a new row for the final dataframe
    rf_row <- process_row(
      rf_fit_pars, intervals[[i]],
      "rf", data_pars_lm_names[[j]], rf_output
    )

    print(paste("Time so far: ", as.numeric(difftime(Sys.time(), start_time, units = "mins")), " minutes"))
    print(rf_row)
    rf_row
  }

  df_rf <- write_results_to_csv(results_rf, prefix)
}

if (all_switch) {
  results_all <- bind_rows(df_optim, df_lm, df_rf) %>%
    arrange(data_pars, basic_pars, data_name)

  write.csv(results_all, "./data/countrywide_results_all.csv")
}