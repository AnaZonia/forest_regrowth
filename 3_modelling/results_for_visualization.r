


# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)

# Source other scripts
source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_normalize_cross_validate.r")
source("3_modelling/2_feature_selection.R")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------ Calculate Permutation Importance ---------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

calculate_permutation_importance <- function(model, data, data_pars) {
    data_pars <- c(unlist(data_pars), "age")
    # Get baseline performance
    baseline_pred <- growth_curve(model$par, data, model$par["lag"])
    baseline_r2 <- calc_r2(data, baseline_pred)

    # Identify categorical variables
    categorical_vars <- list(
        ecoregion = grep("^ecoreg_", data_pars, value = TRUE),
        topography = grep("^topography_", data_pars, value = TRUE)
    )

    # Calculate importance for each variable
    importance <- numeric(length(data_pars))
    names(importance) <- data_pars

    for (i in seq_along(data_pars)) {
        var_name <- data_pars[i]

        # Skip if already permuted as part of a categorical variable
        if (var_name %in% unlist(categorical_vars)) next

        # Create permuted dataset
        data_permuted <- data

        if (var_name %in% unlist(categorical_vars$ecoregion)) {
            # Permute all ecoregion dummy variables together
            ecoregion_vars <- categorical_vars$ecoregion
            permuted_indices <- sample(nrow(data))
            data_permuted[ecoregion_vars] <- data[ecoregion_vars][permuted_indices, ]
        } else if (var_name %in% unlist(categorical_vars$topography)) {
            # Permute all topography dummy variables together
            topography_vars <- categorical_vars$topography
            permuted_indices <- sample(nrow(data))
            data_permuted[topography_vars] <- data[topography_vars][permuted_indices, ]
        } else {
            # Permute single variable
            data_permuted[[var_name]] <- sample(data[[var_name]])
        }

        # Get performance with permuted variable
        permuted_pred <- growth_curve(model$par, data_permuted, model$par["lag"])
        permuted_r2 <- calc_r2(data_permuted, permuted_pred)

        # Importance = decrease in performance
        importance[i] <- baseline_r2 - permuted_r2
    }

    # Create dataframe for plotting
    importance_df <- data.frame(
        variable = names(importance),
        importance = importance,
        importance_pct = 100 * importance / sum(importance)
    )

    # Sort by importance
    importance_df <- importance_df[order(-importance_df$importance), ]

    return(importance_df)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save trained model for the best parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


data <- import_data("grid_1k_amazon_secondary", biome = 1, n_samples = 10000)
norm_data <- normalize_independently(data)
saveRDS(norm_data$train_stats, file = "./0_results/grid_1k_amazon_secondary_train_stats.rds")
norm_data <- norm_data$train_data

for (basic_pars_name in names(basic_pars_options)) {
    basic_pars <- basic_pars_options[[basic_pars_name]]
    data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
    # data_pars <- c("num_fires", "sur_cover", "dist")
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    saveRDS(model, file = paste0("./0_results/amazon_model_", basic_pars_name, ".rds", sep = ""))
    
    model <- readRDS(paste0("./0_results/amazon_model_lag.rds", sep = ""))
    pred_vs_obs <- data.frame(
        age = norm_data$age,
        pred = growth_curve(model$par, norm_data)
    )

    # get the biomass predictions for the ages 1 to 105 (in 35 year steps)
    for (i in c(35, 70)) {
        norm_data_future <- norm_data
        norm_data_future$age <- norm_data_future$age + i
        pred_future <- growth_curve(model$par, norm_data_future)
        df_future <- data.frame(
            age = norm_data_future$age,
            pred = pred_future
        )
        pred_vs_obs <- rbind(pred_vs_obs, df_future)
    }
    
    if (basic_pars_name == "intercept") {
        pred_vs_obs$pred <- round(pred_vs_obs$pred - model$par["B0"])
    } else {
        # add column obs with the age correspondent to that in norm_data
        pred_vs_obs <- cbind(pred_vs_obs, data.frame(
            obs = norm_data$biomass,
            obs_age = round(norm_data$age + model$par["lag"])
        ))
    }

    write.csv(pred_vs_obs, file = paste0("./0_results/pred_vs_obs_amazon_", basic_pars_name, ".csv"), row.names = FALSE)

}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save feature importance results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


for (asymptote in c("nearest_mature", "full_amazon")) {
    # save permutation importance results for asymptote = "quarter_biomass"
    data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 10000, asymptote = asymptote)
    norm_data <- normalize_independently(data)$train_data

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
    r2 <- calc_r2(norm_data, pred)

    importance_results <- calculate_permutation_importance(model, norm_data, data_pars)

    importance_results$importance_scaled = importance_pct * r2 / 100)

    # write.csv(importance_results, file = paste0("./0_results/importance_", asymptote, ".csv"), row.names = FALSE)
}
