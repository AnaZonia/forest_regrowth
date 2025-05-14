# main.R - Main script that runs experiments

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

library(terra)

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
ncore = 4
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
        result <- run_experiment(experiment, name, 1)
        print(result)
        results <- rbind(results, result)
        write.csv(results, file = "results_yearly_no_aic.csv", row.names = FALSE)
    }
}

# select random 20% of pastureland in the Amazon and apply model predictions

# export CSV with the data for all pastureland. get

# train the model based on historical data

model <- run_optim(train_data, pars_init, conditions)

pred_cv <- growth_curve(model$par, new_future_data, lag = model$par["lag"])

# get predictions with ages = age + 30 (for the next 30 years)

points <- vect(df, geom = c("lon", "lat"), crs = "EPSG:4326")

# export the shapefile with the predicted biomass in the future

# later in GEE convert the shapefile to a raster at the resolution of mapbiomas