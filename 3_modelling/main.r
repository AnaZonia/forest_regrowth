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

# Set up parallel processing
set.seed(1)
ncore = 25
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
data <- import_data("grid_1k_amazon_secondary", biome = biome, n_samples = 10000)

sink("./0_results/amazon_experiments_output.txt")

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

sink()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save trained model for the best parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# biome = Amazon
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("unified_fc", biome = biome, n_samples = 10000)

hist(data$nearest_mature)

norm_data <- normalize_independently(data)
# saveRDS(norm_data$train_stats, file = "./0_results/grid_1k_amazon_secondary_train_stats.rds")
norm_data <- norm_data$train_data


# ---- Helper to Fit and Save Model ----
fit_and_save_model <- function(basic_pars_name) {
    basic_pars <- basic_pars_options[[basic_pars_name]]
    # data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
    data_pars <- c("num_fires", "sur_cover")
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    # saveRDS(model, file = paste0("./0_results/amazon_model_", basic_pars_name, ".rds", sep = ""))
    return(model)
}

# ---- Fit Models ----
model_lag <- fit_and_save_model("lag")

pred_lag <- growth_curve(model_lag$par, norm_data, lag = model_lag$par["lag"])
calc_r2(norm_data, pred_lag)






model_intercept <- fit_and_save_model("intercept")
