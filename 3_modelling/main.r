# main.R - Main script that runs experiments

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

data <- import_data("./0_data/unified_all.csv", biome = biome, n_samples = "all")
coords <- data$coords
data <- data$df


norm_data <- normalize_independently(data)$train_data

basic_pars <- basic_pars_options[["lag"]]
data_pars <- c("sur_cover", "ecoreg", "mean_srad", "mean_aet")

pars_init <- find_combination_pars(basic_pars, data_pars, norm_data)

model <- run_optim(norm_data, pars_init, conditions)

future_data <- norm_data %>%
    mutate(age = age + 30)

standard_ages <- norm_data %>%
    mutate(age = 25)


# make percent error map for 2020

pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])

percent_error <- (pred - norm_data$nearest_biomass) / norm_data$nearest_biomass

# pred_future <- growth_curve(model$par, future_data, lag = model$par["lag"])

pred_standard <- growth_curve(model$par, standard_ages, lag = model$par["lag"])

coords$pred_standard <- pred_standard
coords$nearest_biomass <- standard_ages$nearest_biomass
head(coords)
coords$percentage <- coords$pred_standard / coords$nearest_biomass
coords$percent_error <- percent_error

points <- vect(coords, geom = c("lon", "lat"), crs = "EPSG:4326")

# save points as a shapefile
writeVector(points, "predictions.shp", overwrite = TRUE)
