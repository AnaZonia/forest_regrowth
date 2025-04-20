# main.R - Main script that runs experiments

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

# Source other scripts
source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_normalize_cross_validate.r")
source("2_modelling/2_feature_selection_ga.R")
source("2_modelling/2_perm_importance.r")

# Get configuration
n_samples <- 10000

# Set up parallel processing
set.seed(1)
ncore = 4
registerDoParallel(cores = ncore)

# Load data
biome = 1
n_samples = 10000
data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)

# Fit the model on the full data
norm_data <- normalize_independently(data)$train_data

# OPTIM GA ----------------------------------------------------------------------

# ini_par <- data.frame(B0 = 80, k0 = 1, theta = 2)

# for (j in 2:ncore){
#     ini_par[j, ] <- c(
#         ini_par[1, "B0"] * (1.5 * runif(1) + .5),
#         ini_par[1, "k0"] * (1.5 * runif(1) + .5),
#         ini_par[1, "theta"] * (1.5 * runif(1) + 0.5)
#     )
# }

# conditions <- list('pars[["theta"]] > 10', 'pars[["theta"]] < 0', 'pars[["k0"]] < 0')

# optim_ga(ini_par, norm_data)

# --------------------------------------------------------------------------------

pars_init <- find_combination_pars(basic_pars = basic_pars_options[["lag"]], "data_pars" = data_pars_options(colnames(data))[["land_use_landscape_only"]], norm_data)

final_model <- run_optim(norm_data, pars_init, conditions)





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
        mean_r2 = mean(cv_results)
    ))
}

results <- data.frame()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------- Check Importance of parameters included ---------------------------------#
# biome = Amazon
# basic_pars = intercept
# just vary data_pars
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for (name in names(data_pars_options(colnames(data)))) {
    result <- run_experiment("intercept", name, 1)
    print(result)
    results <- rbind(result)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------- Check Model Forms ---------------------------------#
# biome = Amazon
# data_pars = all_mean_clim
# just vary basic_pars
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for (name in names(basic_pars_options)) {
    result <- run_experiment(name, "all_mean_climate", 1)
    print(result)
    results <- rbind(result)
}
