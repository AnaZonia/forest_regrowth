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


growth_curve <- function(pars, data, lag = NULL) {

    # Define parameters that are not expected to change yearly (not prec or si)
    non_clim_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k
    if ("m_base" %in% names(pars)) {
        pars[["B0"]] <- 0
        k <- rep(pars[["k0"]], nrow(data))
        k <- k + rowSums(sapply(non_clim_pars, function(par) {
            pars[[par]] * data[[par]]
        }, simplify = TRUE)) * (data[["age"]] + lag)
    } else {
        # Add yearly-changing climatic parameters to the growth rate k (if included in the parameter set)

        for (clim_par in intersect(climatic_pars, names(pars))) {
            years <- seq(2019, 1985, by = -1)
            clim_columns <- paste0(clim_par, "_", years)
            k <- k + rowSums(sapply(clim_columns, function(col) pars[[clim_par]] * data[[col]]))
        }

        # Define whether age is an explicit or implicit parameter (to multiply the other parameters by)
        implicit_age <- if (!"age" %in% names(pars)) data[["age"]] else rep(1, nrow(data))
        # Define whether the intercept k0 is to be included in the growth rate k
        k <- if ("k0" %in% names(pars)) pars[["k0"]] * implicit_age else rep(0, nrow(data))

        # for (par in non_clim_pars) {
        #     print(par)
        #     print(pars[[par]])
        #     print(head(data[[par]]))
        # }
        
        k <- k + rowSums(sapply(non_clim_pars, function(par) pars[[par]] * data[[par]] * implicit_age))
    }

    # k[which(k < 1e-10)] <- 1e-10 # Constrains k to avoid negative values
    # Constrains k to avoid negative values
    if ("k0" %in% names(pars)) {
        k[which(k < 0)] <- -log(1 - mean(data[["agbd"]]) / mean(data[["nearest_mature_biomass"]]))
    } else {
        k[which(k < 0)] <- 1e-10
    }

    k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)

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

likelihood <- function(pars, data, conditions) {

    if ("m_base" %in% names(pars)) {
        # Calculate log-normal scaled base using re_base and parameters
        scaled_base <- exp(re_base * pars["sd_base"] + pars["m_base"])

        m_results <- matrix(0, nrow = nrow(data), ncol = length(scaled_base))

        growth_curves <- sapply(scaled_base, function(lag) growth_curve(pars, data, lag))
        residuals <- sweep(growth_curves, 1, data$agbd, "-")

        result <- sum(residuals^2)
    } else {
        result <- sum((growth_curve(pars, data) - data$agbd)^2)
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
    if ("sd_base" %in% names(pars)) {
        conditions <- c(conditions, list('pars["sd_base"] < 0', 'pars["m_base"] < 0'))
    }
    if ("k0" %in% names(pars)) {
        conditions <- c(conditions, list('pars["k0"] < 0'))
    }

    model <- optim(pars, likelihood, data = train_data, conditions = conditions)

    if (is.null(test_data)) {
        return(model)
    } else {
        filtered_test_data <- filter_test_data(train_data, test_data)
        if ("m_base" %in% names(pars)){
            pred <- growth_curve(
                model$par, filtered_test_data,
                exp(re_base * model$par["sd_base"] + model$par["m_base"])
            )
        } else {
            pred <- growth_curve(model$par, filtered_test_data)
        }
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
        for (var in problematic_vars) {
            print(table(train_data[[var]]))
        }
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
#   - The R-squared values from each fold are stored in `rsq_list`, and the best model
#

cross_valid <- function(data, run_function, pars_iter, conditions = NULL) {
    rsq_list <- c()
    pars_list <- list()
    indices <- sample(c(1:6), nrow(data), replace = TRUE)

    for (index in 1:5) {
        # Define the test and train sets
        test_data <- data[indices == index, ]
        train_data <- data[indices != index & indices != 6, ]
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
        rsq_list <- append(rsq_list, model_output$rsq)
        pars_list[[index]] <- model_output$model_par
    }
    
    norm_unseen_data <- normalize_independently(data[indices != 6, ])$train_data
    best_pars_list <- pars_list[[which.max(rsq_list)]]
    best_pars <- as.vector(best_pars_list)
    names(best_pars) <- colnames(best_pars_list)
    unseen_data <- data[indices == 6, ]
    final_output <- run_function(norm_unseen_data, best_pars, conditions, unseen_data)

    # Calculate mean and standard deviation of R-squared across folds
    result <- list(
        rsq = mean(rsq_list),
        rsq_sd = sd(rsq_list),
        rsq_unseen = final_output$rsq,
        pars = best_pars_list
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

    # Select numeric columns for normalization, excluding specified ones
    norm_cols <- c(names(train_data)[!grepl(paste0(c(unlist(categorical), "agbd", "nearest_mature_biomass"), collapse = "|"), names(train_data))])

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
    # Remove columns that are entirely NA
    train_data_norm <- train_data_norm %>% select(where(~ sum(is.na(.)) < nrow(train_data_norm)))

    if (is.null(test_data)) {
        return(list(train_data = train_data_norm))
    } else {
        test_data_norm <- normalize(test_data)
        test_data_norm <- test_data_norm %>% select(where(~ sum(is.na(.)) < nrow(test_data_norm)))
        return(list(train_data = train_data_norm, test_data = test_data_norm))
    }
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
    all_possible_pars <- unique(unlist(c("nearest_mature_biomass", "age", non_data_pars, data_pars)))
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
    row$rsq_unseen <- cv_output$rsq_unseen

    if (is.null(basic_pars_names)) {
        row$basic_pars <- NA
    } else {
        row$basic_pars <- basic_pars_names
    }

    # Define the desired order of columns
    desired_column_order <- c(
        "biome_name", "data_name", "data_pars", "basic_pars", "model_type",
        "rsq", "rsq_sd", "rsq_unseen", "nearest_mature_biomass", "age"
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
        # i <- 2
        # j <- 4
        # k <- 1
        # l <- 3
        i <- iterations_optim$interval[iter]
        j <- iterations_optim$data_par[iter]
        k <- iterations_optim$biome[iter]
        l <- iterations_optim$basic_par[iter]

        data <- normalize_independently(dataframes[[i]][[k]])$train_data
        data_pars_iter <- data_pars[[k]][[j]]
        basic_pars_iter <- basic_pars[[l]]

        # Initialize parameter vector with basic parameters and theta
        all_pars_iter <- c(setNames(
            rep(0, length(data_pars_iter)),
            c(data_pars_iter)
        ))

        all_pars_iter[["theta"]] <- 1
        basic_pars_iter <- c(basic_pars_iter, "theta")
        all_pars_iter["B0"] <- mean(data[["agbd"]])

        if ("age" %in% basic_pars_iter) {
            all_pars_iter["age"] <- 0
        }

        if ("k0" %in% basic_pars_iter) {
            all_pars_iter["k0"] <- -log(1 - mean(data[["agbd"]]) / mean(data[["nearest_mature_biomass"]]))
        }

        if ("m_base" %in% basic_pars_iter) {
            all_pars_iter["m_base"] <- 0
            all_pars_iter["sd_base"] <- 1
        }

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Handle categorical variables by grouping dummy variables together
        for (cat_var in categorical) {
            dummy_indices <- grep(cat_var, data_pars_iter)
            if (length(dummy_indices) > 0) {
                data_pars_iter <- c(data_pars_iter[-dummy_indices], cat_var)
            }
        }

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Initialize the best model with basic parameters
        remaining <- 1:length(data_pars_iter)
        taken <- length(remaining) + 1 # out of the range of values such that remaining[-taken] = remaining for the first iteration

        # best model list
        best <- list(AIC = 0)
        val <- 0
        best[["par"]] <- all_pars_iter[names(all_pars_iter) %in% basic_pars_iter]
        val <- length(basic_pars)

        base_row <- all_pars_iter
        base_row[names(all_pars_iter)] <- NA
        base_row <- c(likelihood = 0, base_row)
        
        should_continue <- TRUE
        # Iteratively add parameters and evaluate the model. Keep only AIC improvements.
        for (i in 1:length(data_pars_iter)) {
            if (!should_continue) break
            optim_remaining_pars <- foreach(j = remaining[-taken]) %dopar% {
                # check for categorical variables (to be included as a group)
                if (data_pars_iter[j] %in% categorical) {
                    inipar <- c(best$par, all_pars_iter[grep(data_pars_iter[j], names(all_pars_iter))])
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

            if (best$AIC == 0 | best_model_AIC < best$AIC) {
                best$AIC <- best_model_AIC
                best$par <- iter_df[best_model, names(all_pars_iter)]
                best$par <- Filter(function(x) !is.na(x), best$par)
                taken <- which(sapply(data_pars_iter, function(x) any(grepl(x, names(best$par)))))
            } else {
                print("No improvement. Exiting loop.")
                should_continue <- FALSE
            }
        }

        ideal_par_combination <- append(ideal_par_combination, list(best$par))
        write_rds(ideal_par_combination, paste0("./data/", name_export, "_ideal_par_combination.rds"))
    }

    return(ideal_par_combination)
}

