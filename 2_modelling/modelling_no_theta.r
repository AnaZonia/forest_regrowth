# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#                 Forest Regrowth Model Functions and Utilities
#
#                            Ana Avila - August 2024
#
#     This script defines the core functions used in the forest regrowth
#     modeling process (fit_3_run_model.r)
#
#     Functions included:
#     - growth_curve
#     - likelihood
#     - run_optim
#     - run_lm
#     - filter_test_data
#     - calc_r2
#     - cross_valid
#     - process_row
#     - find_combination_pars
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------- Optimization for Forest Regrowth ---------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   Prepares and executes the optimization process for the forest regrowth model.
#   It applies constraints, performs optimization, and optionally calculates R-squared.
#
# Arguments:
#   train_data : Data frame containing the training dataset.
#   pars       : Named vector of initial parameter values for optimization.
#   conditions : List of parameter constraints expressed as character strings.


run_optim <- function(train_data, pars, conditions) {
    if ("B0" %in% names(pars)) {
        conditions <- c(conditions, list('pars["B0"] < 0'))
    } else {
        conditions <- c(conditions, list('pars["lag"] < 0'))
    }

    model <- optim(pars, likelihood, data = train_data, conditions = conditions)

    return(model)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------- Chapman-Richards Growth Curve ----------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   Calculates the Chapman-Richards growth curve based on provided parameters and data.
#   Incorporates yearly-changing climatic parameters, intercept terms, and forest age.

growth_curve <- function(pars, data, lag = 0) {
    # Define parameters that are not expected to change yearly (not prec or si)
    non_clim_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k
    k <- rep(pars[["k0"]], nrow(data))

    age <- data[["age"]]

    if ("lag" %in% names(pars)) {
        pars[["B0"]] <- 0
        age <- age + lag
        # age <- pmin(age, 1)
    }

    if (length(non_clim_pars) > 0) {
        k <- (k + rowSums(sapply(non_clim_pars, function(par) {
            pars[[par]] * data[[par]]
        }, simplify = TRUE))) * (age)
    }

    # keep original_age column to calculate clim_par

    # Add yearly-changing climatic parameters to the growth rate k (if included in the parameter set)
    for (clim_par in intersect(climatic_pars, names(pars))) {
        for (yrs in 1:max(data[["age"]])) {
            indices <- which(data[["age"]] == yrs)
            # Generate a sequence of years for the current age group
            # Starting from 2019 and going back 'yrs' number of years
            last_year <- max(2019 - yrs - round(lag) + 1, 1985)
            year_seq <- seq(2019, last_year, by = -1)
            clim_columns <- paste0(clim_par, "_", year_seq)
            # as.matrix(t()) is used to ensure that rowSums would work in cases with a single row
            k[indices] <- k[indices] + rowSums(as.matrix(t(sapply(clim_columns, function(col) pars[[clim_par]] * data[[col]][indices]))))
        }
    }

    # Constrains k to avoid negative values
    k[which(k < 1e-10)] <- 1e-10
    k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)

    return(pars[["B0"]] + (data[["nearest_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^theta_val) # ^pars[["theta"]])
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------ Likelihood Function for Optimization -------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   Calculates the likelihood (or sum of squared errors) for the growth curve model,
#   incorporating parameter constraints.
#
# Arguments:
#   pars       : Named vector of parameter values to be evaluated.
#   data       : Dataframe containing predictor variables and forest attributes
#   conditions : List of parameter constraints expressed as character strings.
#
# Returns:
#   Numeric value representing the likelihood or sum of squared errors.
#   Returns -Inf if constraints are violated or if the result is invalid.
#
# External Functions:
#   growth_curve()

likelihood <- function(pars, data, conditions) {
    if ("lag" %in% names(pars)) {
        growth_curves <- growth_curve(pars, data, pars[["lag"]])
        residuals <- growth_curves - data$biomass
        result <- mean(residuals^2)
    } else {
        result <- mean((growth_curve(pars, data) - data$biomass)^2)
    }

    # Check whether any of the parameters is breaking the conditions (e.g. negative values)
    if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
        return(-Inf)
    } else if (is.na(result) || result == 0) {
        return(-Inf)
    } else {
        return(result)
    }
}



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
    exclusion_list <- c(unlist(categorical), unlist(paste0(climatic_pars, "_")), unlist(binary), "biomass", "nearest_biomass", "age")
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
        test_data <- test_data %>% filter(rowSums(test_data[clim_cols]) > 0)
        # test_data <- test_data %>% select(where(~ sum(is.na(.)) < nrow(test_data)))
        return(list(train_data = train_data_norm, test_data = test_data, train_stats = train_stats))
    }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------- Identify Optimal Parameter Combination -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   This function identifies the optimal combination of parameters for a given dataset
#   by iteratively fitting parameter combinations with run_optim and selecting the one that minimizes
#   the Akaike Information Criterion (AIC).
#
# Arguments:
#   iterations       : A dataframe where each row contains the information needed to perform
#                      one iteration of the parameter selection process, including the land use history
#                      interval, data parameter set, basic parameter set, and biome.
#
# Returns:
#   ideal_par_combination : A list where each element contains the best parameter combination
#                           identified for each iteration in the `iterations` dataframe.
#
# Notes:
#   - Categorical variables are handled by grouping their dummy variables together during optimization.
#   - The function writes the results of each iteration to an RDS file for future use.
# External Functions:
#   run_optim()


find_combination_pars <- function(basic_pars, data_pars, data) {
    # data <- train_data
    # Initialize parameter vector with data parameters
    all_pars <- c(setNames(
        rep(0, length(data_pars)),
        c(data_pars)
    ))

    # all_pars[["theta"]] <- 1.5
    all_pars[["k0"]] <- 1

    if ("lag" %in% basic_pars) {
        all_pars[["lag"]] <- 2.5
    } else {
        all_pars[["B0"]] <- mean(data[["biomass"]])
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Handle categorical variables by grouping dummy variables together
    for (cat_var in categorical) {
        dummy_indices <- grep(cat_var, data_pars)
        if (length(dummy_indices) > 0) {
            data_pars <- c(data_pars[-dummy_indices], cat_var)
        }
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Initialize the best model with basic parameters
    remaining <- 1:length(data_pars)
    taken <- length(remaining) + 1 # out of the range of values such that remaining[-taken] = remaining for the first iteration

    # best model list
    best <- list(AIC = 0)
    best[["par"]] <- all_pars[names(all_pars) %in% basic_pars]
    val <- 0

    base_row <- all_pars
    base_row[names(all_pars)] <- NA
    base_row <- c(likelihood = 0, base_row)

    should_continue <- TRUE
    # Iteratively add parameters and evaluate the model. Keep only AIC improvements.
    for (i in 1:length(data_pars)) {
        if (!should_continue) break

        # iter_df <- tibble()
        iter_df <- foreach(j = remaining[-taken]) %dopar% {
            # for (j in remaining[-taken]) {
            # print(j)
            # j = 1
            # check for categorical variables or yearly climatic variables (to be included as a group)
            if (data_pars[j] %in% categorical) {
                inipar <- c(best$par, all_pars[grep(data_pars[j], names(all_pars))])
                # } else if (data_pars[j] %in% climatic) {}
                # inipar <- c(best$par, all_pars[grep(data_pars[j], names(all_pars))])
            } else {
                # as starting point, take the best values from last time
                inipar <- c(best$par, all_pars[data_pars[j]])
            }

            model <- run_optim(data, inipar, conditions)
            iter_row <- base_row
            iter_row[names(inipar)] <- model$par
            iter_row["likelihood"] <- model$value

            # iter_df <- bind_rows(iter_df, iter_row)
            return(iter_row)
        }

        iter_df <- as.data.frame(do.call(rbind, iter_df))

        best_model <- which.min(iter_df$likelihood)
        # best_model_AIC <- 2 * iter_df$likelihood[best_model] + 2 * (i + val + 1)
        best_model_AIC <- iter_df$likelihood[best_model]
        print(best_model_AIC)
        if (best$AIC == 0 | best_model_AIC < best$AIC) {
            best$AIC <- best_model_AIC
            best$par <- iter_df[best_model, names(all_pars)]
            best$par <- Filter(function(x) !is.na(x), best$par)
            taken <- which(sapply(data_pars, function(x) any(grepl(x, names(best$par)))))
            print(paste0("num parameters included: ", i, "parameters taken: ", toString(data_pars[taken])))
        } else {
            not_taken <- data_pars[!data_pars %in% names(best$par)]
            print(paste("No improvement. Exiting loop. Parameters not taken:", toString(not_taken)))
            should_continue <- FALSE
        }
    }
    return(best$par)
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
    indices <- sample(c(1:5), nrow(dataframe), replace = TRUE)
    dataframe$pred_cv <- NA
    dataframe$pred_final <- NA
    r2_list <- numeric(5)

    for (index in 1:5) {
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
        pred_cv <- growth_curve(model$par, test_data)

        # save the predicted values of each iteration of the cross validation.
        dataframe$pred_cv[indices == index] <- pred_cv
        r2 <- calc_r2(dataframe[indices == index, ], pred_cv)
        r2_list[index] <- r2
        print(r2)
    }

    return(r2_list)
}





# calculate_permutation_importance <- function(model, data, data_pars) {
#   # Get baseline performance
#   baseline_pred <- growth_curve(final_model$par, data)
#   baseline_r2 <- calc_r2(data, baseline_pred)

#   # Calculate importance for each variable
#   importance <- numeric(length(data_pars))
#   names(importance) <- data_pars



#   for (i in seq_along(data_pars)) {
#     var_name <- data_pars[i]

#     # Create permuted dataset
#     data_permuted <- data
#     data_permuted[[var_name]] <- sample(data[[var_name]])

#     # pars_init <- find_combination_pars(basic_pars = c("B0", "k0", "theta"), "data_pars" = c("num_fires", "sur_cover"), data_permuted)
#     # model_iter <- run_optim(data_permuted, pars_init, conditions)

#     # Get performance with permuted variable
#     permuted_pred <- growth_curve(model$par, data_permuted)
#     permuted_r2 <- calc_r2(data_permuted, permuted_pred)

#     # Importance = decrease in performance
#     importance[i] <- baseline_r2 - permuted_r2
#   }

#   # Create dataframe for plotting
#   importance_df <- data.frame(
#     variable = names(importance),
#     importance = importance,
#     importance_pct = 100 * importance / sum(importance)
#   )

#   # Sort by importance
#   importance_df <- importance_df[order(-importance_df$importance), ]

#   return(importance_df)
# }



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