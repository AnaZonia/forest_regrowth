
library(foreach)
library(doParallel)
library(tidyverse)
library(mgcv)
library(randomForest)
library(car)

setwd("/home/aavila/forest_regrowth")

# Source external R scripts for data import and function definitions
source("./2_R_scripts/1_import_data.r")
source("./2_R_scripts/2_functions.r")

set.seed(1)
ncores <- 50
registerDoParallel(cores = ncores)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

name_import <- "non_aggregated"
name_export <- name_import

# number of rows to be included in analysis
n_samples <- 10000

# List of climatic parameters that change yearly
climatic_pars <- c("prec", "si")
categorical <- c("ecoreg", "soil")

conditions <- list(
    function(pars) pars["theta"] > 10,
    function(pars) pars["theta"] < 0,
    function(pars) pars["B0"] < 0
)


biomes <- c("amaz", "atla", "both")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define land-use history intervals to import four dataframes
intervals <- list("5yr", "10yr", "15yr", "all")

datafiles <- paste0("./data/", name_import, "_", intervals, ".csv")
dataframes <- lapply(datafiles, import_data, convert_to_dummy = TRUE)
dataframes_lm <- lapply(datafiles, import_data, convert_to_dummy = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# check multicollinearity when adding ecoregion and/or soils
# sometimes you can transform collinear variables - can we control from time since last fire
# that would be of interest in relation with age
pars <- c(
    "age", "nearest_mature", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
    "lulc_sum_40", "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
    "cwd"
)

iterations <- expand.grid(
    interval = seq_along(intervals),
    biome = seq_along(biomes)
)

for (iter in 1:nrow(iterations)) {
    i <- iterations$interval[iter]
    k <- iterations$biome[iter]

    data <- dataframes[[i]][[k]]
    print(biomes[[k]])
    
    initial_pars <- find_combination_pars(pars, data)
    optim_cv_output <- cross_valid(data, run_optim, initial_pars, conditions)
    optim_cv_output$rsq

    lm_cv_output <- cross_valid(data, run_lm, pars)
    lm_cv_output$rsq
    lm_output <- lm(data, pars)
}



