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

test_switch <- TRUE # Set to TRUE for single row or FALSE for all rows

all_switch <- TRUE
# optim_switch <- TRUE
optim_switch <- FALSE
# lm_switch <- TRUE
lm_switch <- FALSE
# rf_switch <- TRUE
rf_switch <- FALSE

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
datafiles <- paste0("./data/", intervals, "_LULC_mat_dist.csv")
dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to create training and testing samples

sample_data <- function(df, size_train, size_test) {
  train_dataframes <- sample_n(df, size = size_train)
  remaining <- df[!(row.names(df) %in% row.names(train_dataframes)), ]
  test_dataframes <- sample_n(remaining, size = size_test)
  list(train_dataframes, test_dataframes)
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
    Filter(function(x) !any(climatic_pars %in% x), data_pars), # remove sets with prec or si
    function(x) c(x, "mature_biomass")
  ),
  list(
    c("age"), c("mature_biomass"), c("age", "mature_biomass")
  )
)

data_pars_lm_names <- append(data_pars_names, list("age", "mat_biomass", "age_mat_biomass"))

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------------- Run Model --------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

start_time <- Sys.time()
print(start_time)

# Run optimization with growth curve
if (optim_switch || all_switch) {
  results_optim <- run_foreach(iterations_optim, "optim", run_optim, conditions)
}

# Run linear model
if (lm_switch || all_switch) {
  results_lm <- run_foreach(iterations_lm, "lm", run_lm)
}

# Run random forest
if (rf_switch || all_switch) {
  results_rf <- run_foreach(iterations_rf, "rf", run_rf)
}

# Combine all results if all_switch is TRUE
if (all_switch) {
  results_all <- bind_rows(results_optim, results_lm, results_rf) %>%
    arrange(data_pars, basic_pars, data_name)
  write.csv(results_all, "./data/all_results.csv", row.names = FALSE)
}
