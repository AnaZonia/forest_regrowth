# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#              Cross-Validation and R-Squared
#
#                   Ana Avila - August 2025
#
#  Evaluates the model performance using 5-fold cross-validation.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------ Calculate R-squared -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Computes the coefficient of determination (R²).
#'
#' @param data A dataframe containing observed response values (`biomass` column).
#' @param pred A numeric vector of predicted values from the model.
#'
#' @return A numeric value: the R-squared indicating goodness-of-fit.

calc_r2 <- function(data, pred) {
    obs_pred <- lm(data$biomass ~ pred)
    residuals <- summary(obs_pred)$residuals
    sum_res_squared <- sum(residuals^2)
    total_sum_squares <- sum((data$biomass - mean(data$biomass))^2)
    r2 <- 1 - (sum_res_squared / total_sum_squares)

    return(r2)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------- Apply min-max scaling ------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        print(var)
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------- Error propagation ------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Considering the standard deviation of the biomass values, here we make predictions with a biomass value extracted from a distribution with mean biomass and the given standard deviation by ESA CCI.
#'
#' @param data A dataframe containing the full dataset to be split.
#' @param basic_pars List of basic parameters to pass to the model.
#' @param data_pars Vector of predictor names to include in the model.
#' @param conditions Additional conditions to pass to optim.
#'
#' @return list with:
#' - r2 - the final calculated r2 with the 
#' - lag_list - list of lag values predicted in each fold of the cross validation
#' - pars - dataframe with the parameters fit in each iteration
#'



error_prop <- function(data, basic_pars, data_pars, conditions, n_iter = 1000) {
    # Forward selection
    fs_idx <- sample(nrow(data), floor(0.8 * nrow(data)), replace = FALSE)
    fs_train <- data[fs_idx, ]
    fs_test <- data[-fs_idx, ]

    fs_norm <- normalize_independently(fs_train, fs_test)
    fs_train_norm <- fs_norm$train_data
    fs_train_stats <- fs_norm$train_stats

    fs_result <- forward_selection(basic_pars, data_pars, fs_train_norm)
    selected_pars <- fs_result[[1]] # parameter structure carried into loop
    r2_progression <- fs_result[[2]] # predictor-by-predictor R² table

    # Bootstrap Monte Carlo (×n_iter)
    # Each iteration gets its own split + normalization.
    pars_list <- vector("list", n_iter)
    r2_vec <- numeric(n_iter)

    for (i in seq_len(n_iter)) {
        train_idx <- sample(nrow(data), floor(0.8 * nrow(data)), replace = FALSE)
        train_data <- data[train_idx, ]
        test_data <- data[-train_idx, ]

        # Normalize using THIS iteration's training min/max
        norm_out <- normalize_independently(train_data, test_data)
        train_norm <- norm_out$train_data
        test_norm <- norm_out$test_data

        # Subsample + Monte Carlo
        sample_idx <- sample(nrow(train_norm), 10000, replace = TRUE)
        data_sampled <- train_norm[sample_idx, ]
        data_sampled$biomass <- rnorm(
            nrow(data_sampled),
            mean = data_sampled$biomass,
            sd   = data_sampled$sd
        )

        # Optim starting from the forward-selected structure
        model <- run_optim(data_sampled, selected_pars, conditions)
        pars_df <- as.data.frame(t(model$par))

        # Evaluate on normalized test set using THIS iteration's parameters
        lag_val <- if ("lag" %in% names(pars_df)) pars_df[["lag"]] else 0
        pred <- growth_curve(pars_df, test_norm, lag = lag_val)
        r2_vec[i] <- calc_r2(test_norm, pred)
        pars_list[[i]] <- pars_df

        if (i %% 100 == 0) {
            message(sprintf("Iter %4d / %d  |  R² = %.3f", i, n_iter, r2_vec[i]))
        }
    }

    list(
        r2             = r2_vec,
        pars           = bind_rows(pars_list),
        r2_progression = r2_progression # forward selection R2 increase by parameter inclusion
    )
}
