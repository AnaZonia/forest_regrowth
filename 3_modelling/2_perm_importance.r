
calculate_permutation_importance <- function(model, data, data_pars) {
    data_pars <- c(unlist(data_pars), "age")
    # Get baseline performance
    baseline_pred <- growth_curve(model$par, data)
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
        permuted_pred <- growth_curve(model$par, data_permuted)
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







# Coefficient based approach
analyze_variable_importance <- function(model, data_pars) {
    # Extract coefficients from model results
    data_pars <- options[["all_no_categorical"]]
    model <- final_model
    coeffs <- model$par[names(model$par) %in% data_pars]

    # Get absolute values to measure magnitude of impact
    importance <- abs(coeffs)

    # Normalize to percentage if desired
    importance_pct <- 100 * importance / sum(importance)

    # Create dataframe for plotting
    importance_df <- data.frame(
        variable = names(importance),
        importance = importance,
        importance_pct = importance_pct
    )

    # Sort by importance
    importance_df <- importance_df[order(-importance_df$importance), ]

    return(importance_df)
}


# Leave one variable out

calculate_loocv_importance <- function(data, basic_pars, data_pars, conditions) {
    # Get baseline model with all variables
    full_data_pars <- data_pars
    pars_init <- find_combination_pars(basic_pars, full_data_pars, data)
    full_model <- run_optim(data, pars_init, conditions)
    full_pred <- growth_curve(full_model$par, data)
    full_r2 <- calc_r2(data, full_pred)

    # Calculate importance for each variable
    importance <- numeric(length(data_pars))
    names(importance) <- data_pars

    for (i in seq_along(data_pars)) {
        var_to_exclude <- data_pars[i]
        reduced_data_pars <- data_pars[data_pars != var_to_exclude]

        # Fit model without this variable
        reduced_pars_init <- find_combination_pars(basic_pars, reduced_data_pars, data)
        reduced_model <- run_optim(data, reduced_pars_init, conditions)
        reduced_pred <- growth_curve(reduced_model$par, data)
        reduced_r2 <- calc_r2(data, reduced_pred)

        # Importance = decrease in performance
        importance[i] <- full_r2 - reduced_r2
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