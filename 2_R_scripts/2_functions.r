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
#     - calc_rsq
#     - cross_valid
#     - process_row
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

growth_curve <- function(pars, data, lag) {
    # Get the names of data_pars to iterate over
    k <- rep(pars["k0"], nrow(data))

    data_pars_names <- names(pars)[!names(pars) %in% basic_pars]

    k <- k + rowSums(sapply(data_pars_names, function(par) {
        pars[[par]] * data[[par]]
    })) * (data["age"] + lag)
    # as mentioned above, best to use means*time...which should also include the "base" growth rate.
    return(data[["nearest_mature_biomass"]] * (1 - exp(-k))^pars[["theta"]])
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
    #   result <- sum((growth_curve(pars, data) - data$agbd)^2)
    scaled_base <- (re_base + pars["m_base"]) * pars["sd_base"]

    m_results <- matrix(0, nrow = nrow(data), ncol = length(scaled_base))
    for (i in 1:length(scaled_base)) # calculate mean across all lags, for each individual
    {
        lag <- scaled_base[i]
        m_results[, i] <- dnorm(growth_curve(pars, lag, data) - data$agbd, sd = par["sd"]) # fills in the likelihood across all individuals (row) for each lag (column)
    }

    # need to get the "probability", take mean of these across all values (i.e., marginalize), and then sum those across all datapoints.
    results <- sum(-log(rowMeans(m_results))) # I think it is negative log likelihoods that we need (assuming we are minimizing), but check that.

    # Check whether any of the parameters is breaking the conditions (e.g. negative values)
    if (any(sapply(conditions, function(cond) cond(pars)))) {
        return(-Inf)
    } else if (is.na(result) || result == 0) {
        return(-Inf)
    } else {
        return(result)
    }
}

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
#   test_data  : Optional. Data frame containing the test dataset for R-squared calculation.
#
# Returns:
#   If test_data is NULL:
#     Returns the full optimization result from optim().
#   If test_data is provided:
#     Returns a list containing:
#       model_par : Optimized parameter values.
#       rsq       : R-squared value of model predictions on filtered test data.
#
# External Functions:
#   likelihood()
#   growth_curve()
#   calc_rsq()
#   filtered_data()


run_optim <- function(train_data, pars, conditions, test_data = NULL) {
    if ("age" %in% names(pars)) {
        conditions <- c(conditions, list(function(pars) pars["age"] < 0))
    }

    model <- optim(pars, likelihood, data = train_data, conditions = conditions)

    if (is.null(test_data)) {
        return(model)
    } else {
        filtered_test_data <- filter_test_data(train_data, test_data)
        pred <- growth_curve(model$par, filtered_test_data)
        rsq <- calc_rsq(filtered_test_data, pred)
        print(paste("R-squared:", rsq))

        return(list(
            model_par = t(model$par),
            rsq = rsq
        ))
    }
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------- Filter Test Data to Match Training Data Range -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   Filters the test dataset to ensure that all rows are within the min-max range of
#   the training dataset for each variable.
#
# Arguments:
#   train_data : The training dataset used to determine valid ranges.
#   test_data  : The test dataset to be filtered.
#
# Returns:
#   filtered_test_data : The subset of test_data where all values are within
#                        the training data's min-max range.

filter_test_data <- function(train_data, test_data) {
    # Identify non-factor columns
    non_factor_columns <- sapply(train_data, is.numeric)
    # Apply min and max only to non-factor columns

    train_min <- sapply(train_data[, non_factor_columns], min)
    train_max <- sapply(train_data[, non_factor_columns], max)

    # Function to check if a row is within the range
    is_within_range <- function(row) {
        all(row >= train_min & row <= train_max)
    }

    # Apply the function to each row of test_data
    filtered_test_data <- test_data[apply(test_data[, non_factor_columns], 1, is_within_range), ]

    return(filtered_test_data)
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
    obs_pred <- lm(data$agbd ~ pred)
    residuals <- summary(obs_pred)$residuals
    sum_res_squared <- sum(residuals^2)
    total_sum_squares <- sum((data$agbd - mean(data$agbd))^2)
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
#   - The R-squared values from each fold are stored in `r_squared_fit`, and the best model
#

cross_valid <- function(data, run_function, pars_iter, conditions = NULL) {
    r_squared_fit <- c()
    pars_fit <- list()
    indices <- sample(c(1:5), nrow(data), replace = TRUE)

    for (index in 1:5) {
        # Define the test and train sets
        test_data <- data[indices == index, ]
        train_data <- data[!indices == index, ]

        # Normalize training and test sets independently, but using training data's min/max for both
        norm_data <- normalize_independently(train_data, test_data)
        train_data <- norm_data$train_data
        test_data <- norm_data$test_data

        # Run the model function on the training set and evaluate on the test set
        if (identical(run_function, run_optim)) {
            model_output <- run_function(train_data, pars_iter, conditions, test_data)
        } else {
            model_output <- run_function(train_data, pars_iter, test_data)
        }

        # Collect R-squared values and model parameters
        r_squared_fit <- append(r_squared_fit, model_output$rsq)
        pars_fit[[index]] <- model_output$model_par
    }

    # Calculate mean and standard deviation of R-squared across folds
    result <- list(
        rsq = mean(r_squared_fit),
        rsq_sd = sd(r_squared_fit),
        pars = pars_fit[[which.max(r_squared_fit)]]
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

normalize_independently <- function(train_data, test_data) {
    # Columns to exclude from normalization
    exclude_cols <- c("soil", "biome", "ecoreg", "last_LU", "protec", "indig", "agbd", "nearest_mature_biomass", "fallow")

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
    test_data_norm <- normalize(test_data)

    # Remove columns that are entirely NA (optional)
    train_data_norm <- train_data_norm %>% select(where(~ sum(is.na(.)) < nrow(train_data_norm)))
    test_data_norm <- test_data_norm %>% select(where(~ sum(is.na(.)) < nrow(test_data_norm)))

    return(list(train_data = train_data_norm, test_data = test_data_norm))
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
