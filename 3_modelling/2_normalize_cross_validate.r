# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------- Calculate R-squared -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   Calculates the R-squared value.
#
# Arguments:
#   data : A dataframe containing the observed values.
#   pred : A vector of predicted values from the model.
#
# Returns:
#   r2  : The R-squared value indicating the goodness of fit.

calc_r2 <- function(data, pred) {
    obs_pred <- lm(data$biomass ~ pred)
    residuals <- summary(obs_pred)$residuals
    sum_res_squared <- sum(residuals^2)
    total_sum_squares <- sum((data$biomass - mean(data$biomass))^2)
    r2 <- 1 - (sum_res_squared / total_sum_squares)

    return(r2)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------- Prepare Dataframes Function --------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Function used in import_data to normalize numeric columns in dataframes.
# Arguments:
#   train_data  : The training dataframe to be used for analysis
#   test_data   : (Optional) The testing dataframe to be normalized using the same
#                 parameters as train_data.
# Returns:
#   A list containing:
#     - train_data : A dataframe with normalized numerical values
#     - test_data  : (If provided) The test dataframe, normalized using train_data statistics

normalize_independently <- function(train_data, test_data = NULL) {

    # Identify numeric columns to normalize (excluding those in exclusion_list)
    exclusion_list <- c(categorical, paste0(climatic_pars, "_"), binary, "biomass", "asymptote", "age")

    norm_cols <- names(train_data)[!grepl(paste0(exclusion_list, collapse = "|"), names(train_data))]
    
    # Compute summary statistics
    train_stats <- data.frame(
        variable = norm_cols,
        min = sapply(norm_cols, function(var) min(train_data[[var]], na.rm = TRUE)),
        max = sapply(norm_cols, function(var) max(train_data[[var]], na.rm = TRUE))
    )
    
    # Apply Min-Max scaling using the precomputed min and max
    for (i in seq_along(norm_cols)) {
        var <- norm_cols[i]
        train_data[[var]] <- (train_data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])

        if (!is.null(test_data)) {
            test_data[[var]] <- (test_data[[var]] - train_stats$min[i]) /
                (train_stats$max[i] - train_stats$min[i])
        }
    }

    
    # Compute normalization statistics for each climatic variable across all years
    for (clim_par in climatic_pars) {
        clim_cols <- names(train_data)[grepl(paste0(clim_par, "_"), names(train_data))]
        if (length(clim_cols) == 0) {
            next
        }
        # Compute summary statistics
        clim_stats <- data.frame(
            variable = paste0(clim_par, "_"),
            min = min(as.matrix(train_data[clim_cols]), na.rm = TRUE),
            max = max(as.matrix(train_data[clim_cols]), na.rm = TRUE)
        )

        train_data[clim_cols] <- (train_data[clim_cols] - clim_stats$min) / (clim_stats$max - clim_stats$min)

        if (!is.null(test_data)) {
            test_data[clim_cols] <- (test_data[clim_cols] - clim_stats$min) / (clim_stats$max - clim_stats$min)
        }
    }

    if (is.null(test_data)) {
        return(list(train_data = train_data, train_stats = train_stats))
    } else {
        # keep in test_data only rows with values greater than zero (those with values in the range of the training data)
        return(list(train_data = train_data, test_data = test_data, train_stats = train_stats))
    }
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------------------- K-Fold Cross-Validation ---------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Function Description:
#   This function performs k-fold cross-validation (with k=5) on the provided data.
#   It splits the data into five folds, uses each fold as a test set once while training
#   on the remaining folds, and then applies the specified `run_function` to train the model
#   and calculate the R-squared values for each fold.
#
# Arguments:
#   data        : The full dataset to be split into training and test sets for cross-validation.
#   run_function: The function used to train the model and generate predictions (either "run_optim", "run_lm")
#   pars_iter   : The parameters to be passed to `run_function` for model training.
#   conditions  : Additional conditions or constraints to be passed to `run_function`
#                 (optional, default is NULL).
#
# Returns:
#   A list with the following elements:
#     - `r2`     : The mean R-squared value across all folds.
#     - `r2_sd`  : The standard deviation of the R-squared values across all folds.
#     - `pars`    : The model parameters corresponding to the fold with the highest R-squared value.
#
# Notes:
#   - The function uses random sampling to assign data points to each of the five folds.
#   - The function assumes that `run_function` takes as input the training data, parameters,
#     conditions, and test data, and returns the model output.
#   - The R-squared values from each fold are stored in `r2_list`, and the best model
#


cross_validate <- function(dataframe, basic_pars, data_pars, conditions) {
    # dataframe <- data
    # basic_pars <- basic_pars_options[["lag"]]
    # data_pars <- c("num_fires", "sur_cover")

    indices <- sample(c(1:5), nrow(dataframe), replace = TRUE)
    dataframe$pred_cv <- NA
    dataframe$pred_final <- NA
    r2_list <- numeric(5)
    print(nrow(dataframe))

    for (index in 1:5) {
        # index <- 1
        # Define the test and train sets
        test_data <- dataframe[indices == index, -grep("pred", names(dataframe))]
        train_data <- dataframe[indices != index, -grep("pred", names(dataframe))]

        # Normalize training and test sets independently, but using training data's min/max for both
        norm_data <- normalize_independently(train_data, test_data)
        train_data <- norm_data$train_data
        test_data <- norm_data$test_data

        # Function to perform direct optimization
        pars_init <- find_combination_pars(basic_pars, data_pars, train_data)

        # Run the model function on the training set and evaluate on the test set
        model <- run_optim(train_data, pars_init, conditions)
        
        pred_cv <- growth_curve(model$par, test_data, lag = model$par["lag"])

        # save the predicted values of each iteration of the cross validation.
        dataframe$pred_cv[indices == index] <- pred_cv
        r2 <- calc_r2(dataframe[indices == index, ], pred_cv)
        r2_list[index] <- r2
        print(r2)
    }

    return(r2_list)
}
