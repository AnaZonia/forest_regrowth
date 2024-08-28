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
#     - filter_test_data
#     - run_optim
#     - run_gam
#     - run_lm
#     - run_rf
#     - process_row
#     - run_foreach
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


growth_curve <- function(pars, data) {

    # Define parameters that are not expected to change yearly (not prec or si)
    non_clim_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars))

    # Define whether age is an explicit or implicit parameter (to multiply the other parameters by)
    implicit_age <- if (!"age" %in% names(pars)) data[["age"]] else rep(1, nrow(data))
    # Define whether the intercept k0 is to be included in the growth rate k
    k <- if ("k0" %in% names(pars)) pars[["k0"]] * implicit_age else rep(0, nrow(data))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k
    if (length(non_clim_pars) > 0) {
        k <- k + rowSums(sapply(non_clim_pars, function(par) pars[[par]] * data[[par]] * implicit_age))
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Add yearly-changing climatic parameters to the growth rate k (if included in the parameter set)
    for (clim_par in intersect(climatic_pars, names(pars))) {
        years <- seq(2019, 1985, by = -1)
        clim_columns <- paste0(clim_par, "_", years)
        k <- k + rowSums(sapply(clim_columns, function(col) pars[[clim_par]] * data[[col]]))
    }

    # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)
    k[which(k > 7)] <- 7 
    # Constrains k to avoid negative values
    if ("k0" %in% names(pars)) {
        k[which(k < 0)] <- -log(1 - mean(data[["agbd"]]) / mean(data[["nearest_mature"]]))
    } else {
        k[which(k < 0)] <- 0
    }

    if (fit_logistic) {
        return(data[["nearest_mature"]] * (1 / (1 + exp(k))))
    } else {
        if ("B0" %in% names(pars)) {
            return(pars[["B0"]] + (data[["nearest_mature"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
        } else {
            return(data[["nearest_mature"]] * (1 - exp(-(pars[["B0_exp"]] + k)))^pars[["theta"]])
        }
    }
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
    result <- sum((growth_curve(pars, data) - data$agbd)^2)

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
        conditions <- c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
    }
    if ("B0" %in% names(pars)) {
        conditions <- c(conditions, list('pars["B0"] < 0'))
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
# ---------------------------- Run Linear Model and Evaluate ----------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   Fits a linear model using the specified parameters, handles rank deficiency,
#   predicts on filtered test data, and calculates R-squared.
#
# Arguments:
#   train_data : Data frame containing the training dataset.
#   pars       : Vector of parameter names to be used in the model.
#   conditions : Placeholder - not used here, but used to prevent loops in cross_valid and run_foreach
#   test_data  : Data frame containing the test dataset.
#
# Returns:
#   list containing:
#     model_par : Named matrix / array of model coefficients (excluding intercept).
#     rsq       : R-squared value of the model predictions on filtered test data.
# External Functions:
#   filtered_data()
#   calc_rsq()

run_lm <- function(train_data, pars, test_data) {
    lm_formula <- as.formula(paste("agbd ~", paste(pars, collapse = " + ")))

    model <- lm(lm_formula, data = train_data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Check for rank deficiency
    aliased_vars <- summary(model)$aliased

    if (any(aliased_vars)) {
        problematic_vars <- names(aliased_vars)[aliased_vars]
        print(paste("Rank-deficient variables:", paste(problematic_vars, collapse = ", ")))
    }

    filtered_test_data <- filter_test_data(train_data, test_data)
    pred <- predict(model, newdata = filtered_test_data)
    rsq <- calc_rsq(filtered_test_data, pred)
    print(paste("R-squared:", rsq))

    return(list(
        model_par = t(summary(model)$coefficients[-1, 1, drop = FALSE]), # -1 to remove (Intercept),
        rsq = rsq
    ))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random Forest
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_rf <- function(train_data, pars, test_data) {
    rf_formula <- as.formula(paste("agbd ~", paste(pars, collapse = " + ")))

    model <- randomForest(rf_formula,
        data = train_data,
        ntree = 100, mtry = 2, importance = TRUE,
        keep.forest = TRUE, oob.score = TRUE, do.trace = 10, parallel = TRUE
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Output R-squared value and model results
    filtered_test_data <- filter_test_data(train_data, test_data)
    pred <- predict(model, newdata = filtered_test_data)
    rsq <- cor(filtered_test_data$agbd, pred)^2
    print(paste("R-squared:", rsq))

    return(list(
        model_par = t(importance(model)[, 1]),
        rsq = rsq
    ))
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
        # Run the model function on the training set and evaluate on the test set
        if (identical(run_function, run_optim)){
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------------------ Process Model Output into DataFrame Row --------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   This function processes the output from cross-validation or model fitting into a
#   structured data frame row. It organizes the model parameters, R-squared values, and
#   other relevant information, ensuring missing parameters are filled with NA, and columns
#   are ordered according to the specified structure.
#
# Arguments:
#   cv_output        : List containing cross-validation output, including coefficients and R-squared values.
#   model_type       : Type of model ("optim", "lm", "rf", "gam").
#   data_name        : Name representing the intervals of land use history included (e.g., "5y", "10y", "15y", "all").
#   data_pars_names  : Names of the fit parameters corresponding to data
#   basic_pars_names : Names of the basic parameters used in the model (only for "optim" models)
#   biome_name       : Name of the biome ("amaz", "atla", or "both")
#
# Returns:
#   A single-row data frame containing the organized model output.


# Define helper functions
process_row <- function(
    cv_output, model_type, data_name, data_pars_names, biome_name, basic_pars_names = NULL) {
    # Initialize a data frame with model parameters (coefficients or variable importance)
    row <- as.data.frame(cv_output$pars)
    # Identify parameters missing from this iteration and add them as NA columns
    all_possible_pars <- unique(unlist(c("nearest_mature", "age", non_data_pars, data_pars)))
    missing_cols <- setdiff(all_possible_pars, names(row))
    row[missing_cols] <- NA

    # Reorder the columns to ensure consistent output structure
    row <- row[, all_possible_pars]

    row$biome_name <- biome_name
    row$data_name <- data_name
    row$model_type <- model_type
    row$data_pars <- data_pars_names
    row$rsq <- cv_output$rsq
    row$rsq_sd <- cv_output$rsq_sd

    if (is.null(basic_pars_names)) {
        row$basic_pars <- NA
    } else {
        row$basic_pars <- basic_pars_names
    }

    # Define the desired order of columns
    desired_column_order <- c(
        "biome_name", "data_name", "data_pars", "basic_pars", "model_type",
        "rsq", "rsq_sd", "nearest_mature", "age"
    )

    row <- row %>%
        select(all_of(desired_column_order), all_of(non_data_pars), everything())

    return(row)
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

find_combination_pars <- function(iterations) {
    ideal_par_combination <- list()

    for (iter in 1:nrow(iterations)) {
        # Extract iteration-specific parameters
        i <- iterations$interval[iter]
        j <- iterations$data_par[iter]
        k <- iterations$biome[iter]

        data <- dataframes[[i]][[k]]
        data_pars_iter <- data_pars[[j]]

        # Initialize parameter vector with basic parameters and theta
        all_pars_iter <- c(setNames(
            rep(0, length(data_pars_iter)),
            c(data_pars_iter)
        ))

        if (!fit_logistic) {
            l <- iterations$basic_par[iter]
            basic_pars_iter <- basic_pars[[l]]

            all_pars_iter[["theta"]] <- 1
            basic_pars_iter <- c(basic_pars_iter, "theta")

            if ("B0" %in% basic_pars_iter) {
                all_pars_iter["B0"] <- mean(data[["agbd"]])
            }

            if ("age" %in% basic_pars_iter) {
                all_pars_iter["age"] <- 0
            }

            if ("k0" %in% basic_pars_iter) {
                if ("B0" %in% basic_pars_iter) {
                    all_pars_iter["k0"] <- -log(1 - mean(data[["agbd"]]) / mean(data[["nearest_mature"]]))
                }
                all_pars_iter["B0_exp"] <- -log(1 - mean(data[["agbd"]]) / mean(data[["nearest_mature"]]))
                all_pars_iter["k0"] <- 0
            }

        }

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle categorical variables by grouping dummy variables together

        categorical <- c("ecoreg", "soil", "last_LU")

        for (cat_var in categorical) {
            dummy_indices <- grep(cat_var, data_pars_iter)
            if (length(dummy_indices) > 0) {
                data_pars_iter <- c(data_pars_iter[-dummy_indices], cat_var)
            }
        }

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Initialize the best model with basic parameters
        remaining <- 1:length(data_pars_iter)
        best <- list(AIC = 0)
        val <- 0
        if (!fit_logistic){
            best[["par"]] <- all_pars_iter[names(all_pars_iter) %in% basic_pars_iter]
            val <- length(basic_pars)
        }
        taken <- length(remaining) + 1 # out of the range of values such that remaining[-taken] = remaining for the first iteration

        base_row <- all_pars_iter
        base_row[names(all_pars_iter)] <- NA
        base_row <- c(likelihood = 0, base_row)

        # Iteratively add parameters and evaluate the model. Keep only AIC improvements.
        for (i in 1:length(data_pars_iter)) {
            optim_remaining_pars <- foreach(j = remaining[-taken]) %dopar% {
                # check for categorical variables (to be included as a group)
                if (data_pars_iter[j] %in% c("last_LU", "ecoreg", "soil")) {
                    all_pars_iter_var <- all_pars_iter[grep(data_pars_iter[j], names(all_pars_iter))]
                    inipar <- c(best$par, all_pars_iter_var)
                } else {
                    # as starting point, take the best values from last time
                    inipar <- c(best$par, all_pars_iter[data_pars_iter[j]])
                }

                model <- run_optim(data, inipar, conditions)

                iter_row <- base_row
                iter_row[names(inipar)] <- model$par
                iter_row["likelihood"] <- model$value
                return(iter_row)
            }

            iter_df <- as.data.frame(do.call(rbind, optim_remaining_pars))
            best_model <- which.min(iter_df$likelihood)
            best_model_AIC <- 2 * iter_df$likelihood[best_model] + 2 * (i + val + 1)

            print(paste0("iteration: ", iter, ", num parameters included: ", i))
            print(best_model_AIC)

            if (best$AIC == 0 | best_model_AIC < best$AIC) {
                best$AIC <- best_model_AIC
                best$par <- iter_df[best_model, names(all_pars_iter)]
                best$par <- Filter(function(x) !is.na(x), best$par)
                taken <- which(sapply(data_pars_iter, function(x) any(grepl(x, names(best$par)))))
            } else {
                print("No improvement. Exiting loop.")
                break
            }
        }

        ideal_par_combination <- append(ideal_par_combination, list(best$par))
        write_rds(ideal_par_combination, paste0("./data/", name_export, "_ideal_par_combination.rds"))
    }

    return(ideal_par_combination)
}
