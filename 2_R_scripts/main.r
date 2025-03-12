# main.R - Main script that runs experiments

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

# Source other scripts
source("2_R_scripts/parameters.r")
source("2_R_scripts/fit_1_data_processing.r")
source("2_R_scripts/fit_2_modelling.r")

# Get configuration
config <- get_config()

# Set up parallel processing
set.seed(1)
registerDoParallel(cores = config$ncores)

# Function to run a single experiment
run_experiment <- function(basic_pars_name, data_pars_name, biome) {
    # Load data
    data <- import_data(
        "./0_data/unified_fc.csv",
        biome = biome,
        n_samples = config$n_samples
    )

    # Get parameters
    basic_pars <- get_basic_pars(basic_pars_name)
    data_pars <- generate_data_pars(data_pars_name, colnames(data))

    # Run cross-validation
    cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

    # Save results
    result_name <- paste(basic_pars_name, data_pars_name, biome, sep = "_")
    save(cv_results, file = paste0("results/", result_name, ".RData"))

    # Return summary
    return(list(
        basic_pars_name = basic_pars_name,
        data_pars_name = data_pars_name,
        biome = biome,
        mean_r2 = cv_results$mean_r2,
        r2_values = cv_results$r2_values
    ))
}


results <- foreach(i = 1:nrow(experiments), .combine = rbind) %dopar% {
    result <- run_experiment(
        experiments$basic_pars_name[i],
        experiments$data_pars_name[i],
        experiments$biome[i]
    )

    data.frame(
        basic_pars_name = result$basic_pars_name,
        data_pars_name = result$data_pars_name,
        biome = result$biome,
        mean_r2 = result$mean_r2
    )
}

# Print summary table
print(results)

# Find the best model
best_model <- results[which.max(results$mean_r2), ]
cat("Best model:\n")
print(best_model)