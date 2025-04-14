# main.R - Main script that runs experiments

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

# Source other scripts
source("2_R_scripts/1_parameters.r")
source("2_R_scripts/1_data_processing.r")
source("2_R_scripts/1_modelling.r")

# Get configuration
n_samples <- 10000

# Set up parallel processing
set.seed(1)
registerDoParallel(cores = 10)

# Load data
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

# write.csv(results, "0_results/model_comparisons.csv")


# get data for plotting
source("2_R_scripts/1_modelling.r")
# Load data
data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)

# Fit the model on the full data
norm_data <- normalize_independently(data)$train_data
pars_init <- find_combination_pars(basic_pars = basic_pars_options[["lag"]], "data_pars" = data_pars_options(colnames(data))[["land_use_landscape_only"]], norm_data)

final_model <- run_optim(norm_data, pars_init, conditions)


data[["pred_lagged"]] <- growth_curve(final_model$par, norm_data, exp(final_model$par[["lag"]]))
data[["pred_nolag"]] <- growth_curve(final_model$par, norm_data) # with no lag, to give the expected values at low ages
data[["pred_future"]] <- growth_curve(final_model$par, norm_data, exp(final_model$par[["lag"]])+35) # with no lag, to give the expected values at low ages

print(mean(data[["pred_lagged"]]))
print(mean(data[["pred_nolag"]]))
print(mean(data[["pred_future"]]))



write.csv(data, paste0(c("0_results/lagged_nolag_unified_data.csv"), collapse = "_"), row.names = FALSE)



# --------------------------------- Plotting ---------------------------------#

# get data for plotting
source("2_R_scripts/1_modelling.r")
source("2_R_scripts/1_data_processing.r")
source("2_R_scripts/1_parameters.r")

biome <- 1
n_samples <- 10000
# Load data

data <- import_data("./0_data/unified_fc_old_biomass.csv", biome = biome, n_samples = n_samples) %>%
    rename(biomass = b1) %>%
    # remove mean_pdsi column
    select(-c(mean_pdsi))

options <- data_pars_options(colnames(data))
options


norm <- normalize_independently(data)
norm_data <- norm$train_data
norm_stats <- norm$train_stats

pars_init <- find_combination_pars(basic_pars = c("B0", "k0", "theta"), "data_pars" = options[["all_no_categorical"]], norm_data)
pars_init


final_model <- run_optim(norm_data, pars_init, conditions)
final_model$par



permutation_importance <- calculate_permutation_importance(final_model, norm_data, options[["all_no_categorical"]])
