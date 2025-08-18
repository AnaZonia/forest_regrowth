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
# Save feature importance results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


for (asymptote in c("full_amazon", "nearest_mature")) { # "quarter_biomass", "ecoreg_biomass",
    # save permutation importance results for asymptote = "quarter_biomass"
    data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000, asymptote = asymptote)
    norm_data <- normalize_independently(data)$train_data

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
    r2 <- calc_r2(norm_data, pred)

    importance_results <- calculate_permutation_importance(model, norm_data, data_pars)

    importance_results$importance_scaled <- importance_results$importance_pct * r2 / 100

    write.csv(importance_results, file = paste0("./0_results/importance_", asymptote, ".csv"), row.names = FALSE)
}

# checking whether the importance results change over different runs
importance_repeats <- data.frame()
for (asymptote in c("full_amazon", "nearest_mature")) {
    for (i in 1:5) {
        data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000, asymptote = asymptote)
        norm_data <- normalize_independently(data)$train_data

        basic_pars <- basic_pars_options[["lag"]]
        data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
        init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
        model <- run_optim(norm_data, init_pars, conditions)
        pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
        r2 <- calc_r2(norm_data, pred)

        importance_results <- calculate_permutation_importance(model, norm_data, data_pars)
        importance_results$importance_scaled <- importance_results$importance_pct * r2 / 100
        importance_wide <- importance_results[, c("variable", "importance_scaled")] %>%
            pivot_wider(names_from = variable, values_from = importance_scaled)
        importance_wide$asymptote <- asymptote
        importance_repeats <- rbind(importance_repeats, importance_wide)
        write.csv(importance_repeats, file = paste0("./0_results/importance_repeats.csv"), row.names = FALSE)
    }
}
