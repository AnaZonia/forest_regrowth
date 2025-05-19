# main.R - Main script
# Runs experiments and saves results for different parameter combinations

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)
library(terra)

# setwd("/home/aavila/Documents/forest_regrowth")

# Source other scripts
source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_normalize_cross_validate.r")
source("3_modelling/2_feature_selection_ga.R")
source("3_modelling/2_perm_importance.r")

# Set up parallel processing
set.seed(1)
ncore = 4
registerDoParallel(cores = ncore)

biome = 1
n_samples = 10000

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check Importance of parameters included
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# biome = Amazon
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

results <- data.frame()

for (name in names(data_pars_options(colnames(data)))) {
    print(data_pars_options(colnames(data))[[name]])
    print("------------------------------------------------")
    for (experiment in names(basic_pars_options)) {
        print(experiment)
        result <- run_experiment(experiment, name, 1)
        print(result)
        results <- rbind(results, result)
        write.csv(results, file = "./0_results/amazon_experiments.csv", row.names = FALSE)
    }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save trained model for the best parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# biome = Amazon
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
