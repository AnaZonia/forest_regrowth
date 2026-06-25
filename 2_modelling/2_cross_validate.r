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
# --------------- 5-Fold Cross-Validation ------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Each fold is used once as a test set while training on
#' the remaining folds. Computes R² values for each fold.
#'
#' @param data A dataframe containing the full dataset to be split.
#' @param basic_pars List of basic parameters to pass to the model.
#' @param data_pars Vector of predictor names to include in the model.
#' @param conditions Additional conditions to pass to optim.
#' @param folds How many times to run the cross-validation (default 5 times)
#'
#' @return list with:
#' - r2_list - list of R² for each one of the five folds
#' - r2_df - dataframe with the R² added per parameter addition (forward selection)
#' - lag_list - list of lag values predicted in each fold of the cross validation
#' - pars - dataframe with the parameters fit in each iteration
#'
#' @details
#' - Each fold is randomly assigned using equal probability.
#' - Training and test sets are normalized independently, with the test set
#'   scaled according to the training set's min/max values.
#' - The model is trained using `run_optim` and evaluated using `calc_r2`.


cross_validate <- function(data, basic_pars, data_pars, conditions, folds = 5) {

    indices <- sample(c(1:folds), nrow(data), replace = TRUE)
    data$pred_cv <- NA
    data$pred_final <- NA
    r2_list <- numeric(folds)
    lag_list <- numeric(folds)
    r2_df <- data.frame()
    pars <- data.frame()

    for (index in 1:folds) {
        # Define the test and train sets
        test_data <- data[indices == index, -grep("pred", names(data))]
        train_data <- data[indices != index, -grep("pred", names(data))]

        # Normalize training and test sets independently, but using training data's min/max for both
        norm_data <- normalize_independently(train_data, test_data)
        train_data <- norm_data$train_data
        test_data <- norm_data$test_data

        # Function to perform direct optimization
        pars_init <- forward_selection(basic_pars, data_pars, train_data)
        # save the R2 increase with each parameter included
        r2_df <- rbind(r2_df, pars_init[[2]])

        # Run the model function on the training set and evaluate on the test set
        model <- run_optim(train_data, pars_init[[1]], conditions)

        pred_cv <- growth_curve(model$par, test_data,
        lag = if ("lag" %in% names(model$par)) model$par["lag"] else 0)

        # save the parameters for each iteration
        pars_df <- as.data.frame(t(model$par))
        print(pars_df)
        pars <- bind_rows(pars, pars_df)

        # save the predicted values of each iteration of the cross validation.
        data$pred_cv[indices == index] <- pred_cv
        r2 <- calc_r2(data[indices == index, ], pred_cv)
        r2_list[index] <- r2
        lag_list[index] <- model$par["lag"]
        print(r2)
        print(model$par["lag"])
    }

    if (nrow(r2_df) > 5) {
        r2_df <- r2_df %>%
            group_by(par) %>%
            summarise(
                mean_r2_diff = mean(r2_diff, na.rm = TRUE),
                sd_r2_diff = sd(r2_diff, na.rm = TRUE),
                n = n()
            )

        r2_df <- r2_df[order(r2_df$mean_r2_diff, decreasing = FALSE), ]
    }

    pars <- pars[, !grepl("lag|k0", names(pars))]

    return(list(r2_list, r2_df, lag_list, pars))
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

error_prop <- function(data, basic_pars, data_pars, conditions) {

    # Define the test and train sets
    ind <- sample(c(TRUE, FALSE), nrow(data), replace = TRUE, prob = c(0.8, 0.2))
    test_data <- data[!ind, ]
    train_data <- data[ind, ]

    # Normalize training and test sets independently, but using training data's min/max for both
    norm_data <- normalize_independently(train_data, test_data)
    train_stats <- norm_data$train_stats
    train_data <- norm_data$train_data
    test_data <- norm_data$test_data

    # Function to perform forward selection
    pars_init <- forward_selection(basic_pars, data_pars, train_data)
    r2_df <- pars_init[[2]] # save the r2 increase by adding each predictor

    pars <- data.frame()

    for (index in 1:10) {

        # Define the test and train sets
        ind <- sample(c(TRUE, FALSE), nrow(data), replace = TRUE, prob = c(0.8, 0.2))
        test_data <- data[!ind, ]
        train_data <- data[ind, ]

        norm_data <- normalize_independently(train_data, test_data)
        train_stats <- norm_data$train_stats
        train_data <- norm_data$train_data
        test_data <- norm_data$test_data


        # randomly select 10,000 rows
        data_sampled <- train_data[sample(nrow(train_data), 10000), ]
        data_sampled$biomass <- rnorm(nrow(data_sampled), data_sampled$biomass, data_sampled$sd)

        # Run the model function on the training set and evaluate on the test set
        model <- run_optim(data_sampled, pars_init[[1]], conditions)

        # save the parameters for each iteration
        pars_df <- as.data.frame(t(model$par))
        print(pars_df)
        pars <- bind_rows(pars, pars_df)

        pred <- growth_curve(pars, test_data,
            lag = if ("lag" %in% names(pars)) final_pars["lag"] else 0
        )

        r2 <- calc_r2(test_data, pred)
    }

    return(list(r2, pars, train_stats))
}





error_prop <- function(data, basic_pars, data_pars, conditions, n_iter = 1000) {
    # ── Phase 1: Forward Selection (runs once on a stable split) ─────────────
    # Keep this split separate and fixed — it only determines model structure,
    # not parameter distributions. Don't reuse it in the bootstrap loop.
    fs_idx <- sample(nrow(data), floor(0.8 * nrow(data)), replace = FALSE)
    fs_train <- data[fs_idx, ]
    fs_test <- data[-fs_idx, ]

    fs_norm <- normalize_independently(fs_train, fs_test)
    fs_train_norm <- fs_norm$train_data

    fs_result <- forward_selection(basic_pars, data_pars, fs_train_norm)
    selected_pars <- fs_result[[1]] # parameter structure carried into loop
    r2_progression <- fs_result[[2]] # predictor-by-predictor R² table

    # ── Phase 2: Bootstrap Monte Carlo (×n_iter) ─────────────────────────────
    # Each iteration gets its own split + normalization. Normalization must be
    # fit on that iteration's training data only, then applied to its test data.
    pars_list <- vector("list", n_iter)
    r2_vec <- numeric(n_iter)
    stats_list <- vector("list", n_iter) # store all norm stats, not just last

    for (i in seq_len(n_iter)) {
        # 1. Bootstrap resample — deterministic size avoids degenerate small splits
        train_idx <- sample(nrow(data), floor(0.8 * nrow(data)), replace = FALSE)
        train_data <- data[train_idx, ]
        test_data <- data[-train_idx, ]

        # 2. Normalize using THIS iteration's training min/max
        norm_out <- normalize_independently(train_data, test_data)
        train_norm <- norm_out$train_data
        test_norm <- norm_out$test_data
        stats_list[[i]] <- norm_out$train_stats

        # 3. Subsample + Monte Carlo noise injection
        #    replace = TRUE allows proper bootstrap resampling of the 10k draw
        sample_idx <- sample(nrow(train_norm), 10000, replace = TRUE)
        data_sampled <- train_norm[sample_idx, ]
        data_sampled$biomass <- rnorm(
            nrow(data_sampled),
            mean = data_sampled$biomass,
            sd   = data_sampled$sd
        )

        # 4. Optimize — always starts from the forward-selected structure
        model <- run_optim(data_sampled, selected_pars, conditions)
        pars_df <- as.data.frame(t(model$par))

        # 5. Evaluate on normalized test set using THIS iteration's parameters
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
        r2_progression = r2_progression, # forward selection diagnostics
        train_stats    = stats_list # all iterations, not just last
    )
}