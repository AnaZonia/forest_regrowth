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
#     - calc_rsq
#     - cross_valid
#     - normalize_independently
#     - find_combination_pars
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------- Chapman-Richards Growth Curve ----------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   Calculates the Chapman-Richards growth curve based on provided parameters and data.
#   Incorporates yearly-changing climatic parameters, intercept terms, and forest age.
#
# Arguments:
#   pars : Named vector of parameters for the growth curve.
#   data : Dataframe containing predictor variables and forest attributes
#
# Returns:
#   Vector of predicted aboveground biomass density (AGBD) values.
#
# Notes:
#   - Supports both logistic and exponential growth models.
#   - Handles forest age as either an explicit parameter (part of "pars")
#       or as an implicit parameter (multiplying all non-yearly predictors by age)
#   - The intercept term can be defined either as B0 or B0_exp (out or in of the exponential)
#   - Incorporates growth rate intercept term k0 if provided
#   - Incorporates yearly-changing climatic parameters if provided.

growth_curve_lag_k <- function(pars, data, lag) {
    k <- rep(pars[["k"]], nrow(data))

    # Get the names of data_pars to iterate over
    data_pars_names <- names(pars)[!names(pars) %in% basic_pars]

    if (length(data_pars_names) > 0) {
        k <- k + rowSums(sapply(data_pars_names, function(par) {
            pars[[par]] * data[[par]]
        }, simplify = TRUE)) * (data[["age"]] + lag)
    }

    k[k <= 0] <- 1e-10

    return(data[["nearest_mature_biomass"]] * (1 - exp(-k))^pars[["theta"]])
}

# growth_curve_lag_A <- function(pars, data, lag) {
#     A <- rep(pars["A"], nrow(data))

#     # Get the names of data_pars to iterate over
#     data_pars_names <- names(pars)[!names(pars) %in% basic_pars]

#     if (length(data_pars_names) > 0) {
#         A <- A + rowSums(sapply(data_pars_names, function(par) {
#             pars[[par]] * data[[par]]
#         }, simplify = TRUE))
#     }

#     return(A * (1 - exp(-pars[["k"]] * (data[["age"]] + lag)))^pars[["theta"]])
# }

# growth_curve_lag_k_A <- function(pars, data, lag) {
#     k <- rep(pars["k"], nrow(data))

#     # Get the names of data_pars to iterate over
#     data_pars_names <- names(pars)[!names(pars) %in% basic_pars]

#     if (length(data_pars_names) > 0) {
#         k <- k + rowSums(sapply(data_pars_names, function(par) {
#             pars[[par]] * data[[par]]
#         }, simplify = TRUE)) * (data[["age"]] + lag)
#     }

#     return(data[["nearest_mature_biomass"]] * (1 - exp(-k))^pars[["theta"]])
# }


growth_curve_B0_theta <- function(pars, data) {
    k <- rep(pars[["k"]], nrow(data))

    # Get the names of data_pars to iterate over
    data_pars_names <- names(pars)[!names(pars) %in% basic_pars]

    if (length(data_pars_names) > 0) {
        k <- k + rowSums(sapply(data_pars_names, function(par) {
            pars[[par]] * data[[par]]
        }, simplify = TRUE)) # * (data[["age"]])
    }

    k[which(k > 7)] <- 7 #
    k[which(k < 0)] <- 1e-10

    return(pars[["B0"]] + (data[["nearest_mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
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


likelihood_lag <- function(pars, data, conditions, growth_curve_func, NLL) {
    # Calculate log-normal scaled base using re_base and parameters
    scaled_base <- exp((re_base + pars["m_base"]) * pars["sd_base"])

    m_results <- matrix(0, nrow = nrow(data), ncol = length(scaled_base))

    growth_curves <- sapply(scaled_base, function(lag) growth_curve_func(pars, data, lag))
    residuals <- sweep(growth_curves, 1, data$agbd, "-")

    if (NLL) {
        m_results <- dnorm(residuals, sd = pars["sd"])
        mean_results <- rowMeans(m_results)
        mean_results[mean_results == 0] <- 1e-10
        # Calculate log-likelihood
        result <- sum(-log(mean_results))
    } else {
        result <- sum(residuals^2)
    }

    # Check parameter constraints
    if (any(sapply(conditions, function(cond) cond(pars)))) {
        return(-Inf)
    } else if (is.na(result) || result == 0) {
        return(-Inf)
    } else {
        return(result)
    }
}

likelihood_B0_theta <- function(pars, data, conditions, growth_curve_func, NLL) {
    # Use NLS if required, otherwise calculate log-likelihood
    residuals <- growth_curve_func(pars, data) - data$agbd

    if (NLL) {
        # Calculate the residuals and the likelihood
        m_results <- dnorm(residuals, sd = pars["sd"])

        # Prevent log(0) by adding a small constant
        mean_results <- mean(m_results)
        mean_results[mean_results == 0] <- 1e-10

        result <- sum(-log(mean_results)) # + 0.01 * log(pars["sd"])

        # if (is.na(result) || is.nan(result)) {
        #     cat("Result is NA or NaN. Debugging information:\n")
        #     cat("pars['sd'] =", pars["sd"], "\n")
        #     cat("mean_results =", mean_results, "\n")
        # } else if (result < 0) {
        #     cat("Result is negative. pars['sd'] =", pars["sd"], "\n")
        # }
        # print(result)

    } else {
        result <- sum(residuals^2)
    }

    # Check parameter constraints
    if (any(sapply(conditions, function(cond) cond(pars)))) {
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
#   rsq  : The R-squared value indicating the goodness of fit.

calc_rsq <- function(data, pred) {
    # Identify indices of infinite values in pred
    infinite_indices <- is.infinite(pred)
    # cat(sprintf("Number of -Inf values output by optim (removed): %d\n", sum(infinite_indices)))
    # Create a logical mask for finite values
    finite_mask <- !infinite_indices
    # Subset pred and data$agbd using the finite mask
    pred_cleaned <- pred[finite_mask]
    agbd_cleaned <- data$agbd[finite_mask]

    obs_pred <- lm(agbd_cleaned ~ pred_cleaned)
    residuals <- summary(obs_pred)$residuals
    sum_res_squared <- sum(residuals^2)
    total_sum_squares <- sum((agbd_cleaned - mean(agbd_cleaned))^2)
    rsq <- 1 - (sum_res_squared / total_sum_squares)

    return(rsq)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------------------- K-Fold Cross-Validation ---------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
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
#     - `rsq`     : The mean R-squared value across all folds.
#     - `rsq_sd`  : The standard deviation of the R-squared values across all folds.
#     - `pars`    : The model parameters corresponding to the fold with the highest R-squared value.
#
# Notes:
#   - The function uses random sampling to assign data points to each of the five folds.
#   - The function assumes that `run_function` takes as input the training data, parameters,
#     conditions, and test data, and returns the model output.
#   - The R-squared values from each fold are stored in `r2_list`, and the best model
#

cross_valid <- function(data, pars, conditions, likelihood_func, growth_curve_func, NLL) {
    r2_list <- c()
    pars_list <- list()
    indices <- sample(c(1:5), nrow(data), replace = TRUE)

    for (index in 1:5) {

        # Define the test and train sets
        test_data <- data[indices == index, ]
        train_data <- data[!indices == index, ]

        # Normalize training and test sets independently, but using training data's min/max for both
        norm_data <- normalize_independently(train_data, test_data)
        train_data <- norm_data$train_data
        test_data <- norm_data$test_data

        model <- optim(pars, function(pars) likelihood_func(pars, train_data, conditions, growth_curve_func, NLL))

        if (identical(likelihood_func, likelihood_B0_theta)) {
            pred <- growth_curve_func(model$par, test_data)
        } else {
            pred <- unname(growth_curve_func(
                model$par, test_data,
                exp((re_base + model$par["m_base"]) * model$par["sd_base"])
            ))
        }

        rsq <- calc_rsq(test_data, pred)

        # Collect R-squared values and model parameters
        r2_list <- append(r2_list, rsq)
        pars_list[[index]] <- model$par
    }

    norm_unseen_data <- normalize_independently(unseen_data)$train_data
    best_pars <- pars_list[[which.max(r2_list)]]
    if (identical(likelihood_func, likelihood_B0_theta)) {
        pred <- growth_curve_func(best_pars, norm_unseen_data)
    } else {
        pred <- unname(growth_curve_func(
            best_pars, norm_unseen_data,
            exp((re_base + best_pars["m_base"]) * best_pars["sd_base"])
        ))
    }
    r2_unseen <- calc_rsq(norm_unseen_data, pred)

    # Calculate mean and standard deviation of R-squared across folds
    result <- list(
        r2_mean = mean(r2_list),
        r2_sd = sd(r2_list),
        r2_unseen = r2_unseen,
        best_pars = best_pars
    )

    return(result)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------- Prepare Dataframes Function --------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Function used in import_data to normalize numeric columns in dataframes.
# Arguments:
#   data             : The dataframe to be used for analysis
# Returns:
#   data             : A dataframe with normalized numerical values

normalize_independently <- function(train_data, test_data = NULL) {
    # Columns to exclude from normalization
    exclude_cols <- c("nearest_mature_biomass", "agbd")

    # Select numeric columns for normalization, excluding specified ones
    norm_cols <- setdiff(names(train_data)[sapply(train_data, is.numeric)], exclude_cols)

    # Compute min and max for normalization based on training data
    train_min_max <- train_data %>%
        summarise(across(all_of(norm_cols), list(min = ~ min(., na.rm = TRUE), max = ~ max(., na.rm = TRUE))))

    # Normalize training and test data using training min/max values
    normalize <- function(data) {
        data %>%
            mutate(across(
                all_of(norm_cols),
                ~ (. - train_min_max[[paste0(cur_column(), "_min")]]) /
                    (train_min_max[[paste0(cur_column(), "_max")]] - train_min_max[[paste0(cur_column(), "_min")]])
            ))
    }

    train_data_norm <- normalize(train_data)
    # Remove columns that are entirely NA (optional)
    # train_data_norm <- train_data_norm %>% select(where(~ sum(is.na(.)) < nrow(train_data_norm)))

    if (is.null(test_data)) {
        return(list(train_data = train_data_norm))
    } else {
        test_data_norm <- normalize(test_data)
        # test_data_norm <- test_data_norm %>% select(where(~ sum(is.na(.)) < nrow(test_data_norm)))
        return(list(train_data = train_data_norm, test_data = test_data_norm))
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

find_combination_pars <- function(pars, data) {
    # Initialize parameter vector with basic parameters and theta
    all_pars_iter <- c(setNames(
        rep(0, length(pars)),
        c(pars)
    ))

    all_pars_iter[["theta"]] <- 1
    all_pars_iter["B0"] <- mean(data[["agbd"]])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Initialize the best model with basic parameters
    remaining <- 1:length(pars)
    taken <- length(remaining) + 1 # out of the range of values such that remaining[-taken] = remaining for the first iteration

    # best model list
    best <- list(AIC = 0)

    basic_pars <- c("theta", "B0")
    # beginning with the essential parameters
    best[["par"]] <- all_pars_iter[names(all_pars_iter) %in% basic_pars]
    val <- length(basic_pars)

    base_row <- all_pars_iter
    base_row[names(all_pars_iter)] <- NA
    base_row <- c(likelihood = 0, base_row)

    # Iteratively add parameters and evaluate the model. Keep only AIC improvements.
    for (i in 1:length(pars)) {
        iter_df <- data.frame()
        optim_remaining_pars <- for (j in remaining[-taken]) {
            # as starting point, take the best values from last time
            print(j)
            inipar <- c(best$par, all_pars_iter[pars[j]])

            model <- run_optim(data, inipar, conditions)

            iter_row <- base_row
            iter_row[names(inipar)] <- model$par
            iter_row["likelihood"] <- model$value
            iter_row <- as.data.frame(t(iter_row))
            iter_df <- rbind(iter_df, iter_row)
        }

        best_model <- which.min(iter_df$likelihood)
        best_model_AIC <- 2 * iter_df$likelihood[best_model] + 2 * (i + val + 1)

        if (best$AIC == 0 | best_model_AIC < best$AIC) {
            best$AIC <- best_model_AIC
            best$par <- iter_df[best_model, names(all_pars_iter)]
            best$par <- Filter(function(x) !is.na(x), best$par)
            taken <- which(sapply(pars, function(x) any(grepl(x, names(best$par)))))
        } else {
            print("No improvement. Exiting loop.")
            break
        }
    }

    return(best$par)
}
