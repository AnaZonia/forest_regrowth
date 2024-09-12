
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

# for (i in 1:length(dataframes_lm)){
#     write.csv(dataframes_lm[[i]][[1]], paste0("processed_amaz_", intervals[[i]], ".csv"))
#     write.csv(dataframes_lm[[i]][[2]], paste0("processed_atla_", intervals[[i]], ".csv"))
# }

pars <- c(
    "age", "nearest_mature", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
    "lulc_sum_40", "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
    "cwd"
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------ Test Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- dataframes[[4]][[2]]

source("./2_R_scripts/2_functions.r")
initial_pars <- find_combination_pars(pars, data)
initial_pars
optim_cv_output <- cross_valid(data, run_optim, initial_pars, conditions)
optim_cv_output$rsq

lm_cv_output <- cross_valid(data, run_lm, c("age", "nearest_mature"))
lm_cv_output$rsq
initial_pars
