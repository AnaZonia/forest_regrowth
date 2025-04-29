# main.R - Main script that runs experiments

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

setwd("/home/aavila/Documents/forest_regrowth")

# Source other scripts
source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_normalize_cross_validate.r")
source("2_modelling/2_feature_selection_ga.R")
source("2_modelling/2_perm_importance.r")

# Get configuration

# Set up parallel processing
set.seed(1)
ncore = 25
registerDoParallel(cores = ncore)

# Load data
biome = 1
n_samples = 10000

data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)

# Function to run a single experiment
run_experiment <- function(basic_pars_name, data_pars_name, biome) {

    # Get parameters
    basic_pars <- basic_pars_options[[basic_pars_name]]
    data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

    # Run cross-validation
    cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

    # Return summary
    return(data.frame(
        basic_pars_name = basic_pars_name,
        data_pars_name = data_pars_name,
        biome = biome,
        mean_r2 = mean(cv_results),
        sd_r2 = sd(cv_results)
    ))
}



# Function to run a single experiment
run_experiment_no_cv <- function(basic_pars_name, data_pars_name, biome) {

    # Get parameters
    basic_pars <- basic_pars_options[[basic_pars_name]]
    data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

    train_data <- normalize_independently(data)$train_data

    pars_init <- find_combination_pars(basic_pars, data_pars, train_data)

    model <- run_optim(train_data, pars_init, conditions)
    pred_cv <- growth_curve(model$par, train_data, lag = model$par["lag"])

    r2 <- calc_r2(data, pred_cv)

    # Return summary
    return(data.frame(
        basic_pars_name = basic_pars_name,
        data_pars_name = data_pars_name,
        r2 =r2
    ))
}


results <- data.frame()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------------- Check Importance of parameters included
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# biome = Amazon
# basic_pars = intercept
# just vary data_pars
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for (name in names(data_pars_options(colnames(data)))) {
    print(data_pars_options(colnames(data))[[name]])
    print("------------------------------------------------")
    for (experiment in names(basic_pars_options)) {
        print(experiment)
        result <- run_experiment_no_cv(experiment, name, 1)
        print(result)
        results <- rbind(results, result)
        write.csv(results, file = "results_no_cv.csv", row.names = FALSE)
    }
}


