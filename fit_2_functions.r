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
#   test_data  : Optional. Data frame containing the test dataset for R-squared calculation.
#
# Returns:
#   If test_data is NULL:
#     Returns the full optimization result from optim().
#   If test_data is provided:
#     Returns a list containing:
#       model_par : Optimized parameter values.
#       r2       : R-squared value of model predictions on filtered test data.
#
# External Functions:
#   likelihood()
#   growth_curve()
#   calc_r2()
#   filtered_data()


run_optim <- function(train_data, pars, conditions) {
    if ("age" %in% names(pars)) {
        conditions <- c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
    }
    if ("B0" %in% names(pars)) {
        conditions <- c(conditions, list('pars["B0"] < 0'))
    }
    if ("sd_base" %in% names(pars)) {
        conditions <- c(conditions, list('pars["sd_base"] < 0', 'pars["m_base"] < 0'))
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
#
# Arguments:
#   pars : Named vector of parameters for the growth curve.
#   data : Dataframe containing predictor variables and forest attributes
#
# Returns:
#   Vector of predicted aboveground biomass density (biomass) values.
#
# Notes:
#   - Supports both logistic and exponential growth models.
#   - Handles forest age as either an explicit parameter (part of "pars")
#       or as an implicit parameter (multiplying all non-yearly predictors by age)
#   - The intercept term can be defined either as B0 or B0_exp (out or in of the exponential)
#   - Incorporates growth rate intercept term k0 if provided
#   - Incorporates yearly-changing climatic parameters if provided.


growth_curve <- function(pars, data, lag = 0) {
    pars <- inipar
    lag = 0
    # Define parameters that are not expected to change yearly (not prec or si)
    non_clim_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars))
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k
    if ("m_base" %in% names(pars)) {
        k <- rep(pars[["k0"]], nrow(data))
        pars[["B0"]] <- 0
        if (length(non_clim_pars) > 0) {
            k <- k + rowSums(sapply(non_clim_pars, function(par) {
                pars[[par]] * data[[par]]
            }, simplify = TRUE)) * (data[["age"]] + lag)
        }
    } else {
        # Define whether age is an explicit or implicit parameter (to multiply the other parameters by)
        implicit_age <- if (!"age" %in% names(pars)) data[["age"]] else rep(1, nrow(data))
        # Define whether the intercept k0 is to be included in the growth rate k
        k <- pars[["k0"]] * implicit_age
        if (length(non_clim_pars) > 0) {
            k <- k + rowSums(sapply(non_clim_pars, function(par) pars[[par]] * data[[par]])) * implicit_age
        }
    }

    # Add yearly-changing climatic parameters to the growth rate k (if included in the parameter set)
    for (clim_par in intersect(climatic_pars, names(pars))) {
        for (yrs in 1:max(data[["age"]])) {
            indices <- which(data[["age"]] == yrs)
            # Generate a sequence of years for the current age group
            # Starting from 2019 and going back 'yrs' number of years
            last_year <- max(2019 - yrs - round(lag) + 1, 1985)
            year_seq <- seq(2019, last_year, by = -1)
            clim_columns <- paste0(clim_par, "_", year_seq)
            k[indices] <- k[indices] + rowSums(sapply(clim_columns, function(col) pars[[clim_par]] * data[indices, col]))
        }
    }
    
    # Constrains k to avoid negative values
    k[which(k < 1e-10)] <- 1e-10
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
        residuals <- sweep(growth_curves, 1, data$biomass, "-")
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

cross_valid <- function(data, pars_iter, conditions = NULL) {

    indices <- sample(c(1:5), nrow(data), replace = TRUE)
    data$pred_cv <- NA
    data$pred_final <- NA
    r2_list <- numeric(5)

    predict_growth <- function(model, data) {
        if ("m_base" %in% names(pars_iter)) {
            growth_curve(
                model$par, data,
                exp(re_base * model$par["sd_base"] + model$par["m_base"])
            )

        } else {
            growth_curve(model$par, data)
        }
    }

    if (is.null(conditions)) { # run lm
        # Filter out elements containing any categorical pattern and keep track of matched categories
        pars_iter_filtered <- names(pars_iter)[!sapply(names(pars_iter), function(name) {
            any(sapply(categorical, function(cat) grepl(cat, name)))
        })]
        pars_iter_filtered <- pars_iter_filtered[!pars_iter_filtered %in% non_data_pars]
        # Identify which categories were found in pars_iter and append only once per matched category
        matched_categories <- categorical[sapply(categorical, function(cat) any(grepl(cat, names(pars_iter))))]
        pars_iter_filtered <- c(pars_iter_filtered, matched_categories, "nearest_mature_biomass")
        lm_formula <- as.formula(paste("biomass ~", paste(pars_iter_filtered, collapse = " + ")))
    }

    for (index in 1:5) {
        # Define the test and train sets
        test_data <- data[indices == index, -grep("pred", names(data))]
        train_data <- data[indices != index, -grep("pred", names(data))]
        # Normalize training and test sets independently, but using training data's min/max for both
        norm_data <- normalize_independently(pars_iter, train_data, test_data)
        train_data <- norm_data$train_data
        test_data <- norm_data$test_data

        if (is.null(conditions)) {  # run lm
            model <- lm(lm_formula, data = train_data)
            pred_cv <- predict(model, test_data)
        } else { # run optim
            # Run the model function on the training set and evaluate on the test set
            model <- run_optim(train_data, pars_iter, conditions)
            pred_cv <- predict_growth(model, test_data)
        }

        # save the predicted values of each iteration of the cross validation.
        data$pred_cv[indices == index] <- pred_cv
        r2 <- calc_r2(data[indices == index, ], pred_cv)
        print(r2)
        r2_list[index] <- r2
    }
    # Fit the model on the full data
    norm_data <- normalize_independently(pars_iter, data[, -grep("pred", names(data))])$train_data

    if (is.null(conditions)) { # run lm
        final_model <- lm(lm_formula, data = norm_data)
        pred_final <- predict(final_model, norm_data)
        data$pred_final <- pred_final
        pars <- coef(final_model)[-1] # remove intercept
        names(pars) <- gsub("([a-zA-Z])([0-9])", "\\1_\\2", names(pars))
    } else {
        final_model <- run_optim(norm_data, pars_iter, conditions)
        pred_final <- predict_growth(final_model, norm_data)
        data$pred_final <- growth_curve(final_model$par, norm_data)
        pars <- as.vector(t(final_model$par))
        names(pars) <- names(final_model$par)
    }

    r2_final <- calc_r2(norm_data, pred_final)

    # Calculate mean and standard deviation of R-squared across folds
    result <- list(
        r2_mean = mean(r2_list, na.rm = TRUE),
        r2_sd = sd(r2_list, na.rm = TRUE),
        r2_final = r2_final,
        pars = pars,
        pred = data$pred_final
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

normalize_independently <- function(pars, train_data, test_data = NULL) {
    # Select numeric columns for normalization, excluding specified ones
    exclusion_list <- c(unlist(categorical), "biomass", "nearest_mature_biomass", "age", climatic_pars)
    # if k is multiplied by the age column, don't normalize age
    exclusion_list <- c(exclusion_list, if (!"age" %in% names(pars)) "age")
    norm_cols <- c(names(train_data)[!grepl(paste0(exclusion_list, collapse = "|"), names(train_data))])

    # Compute mean and standard deviation for normalization based on training data
    train_mean_sd <- train_data %>%
        summarise(across(all_of(norm_cols), list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE))))

    # Compute global mean and standard deviation for each climatic parameter across all years
    climatic_mean_sd <- lapply(climatic_pars, function(param) {
        all_years_values <- train_data %>%
            select(matches(paste0("^", param, "_\\d{4}$"))) %>%
            unlist(use.names = FALSE) # Flatten to a single vector
        list(mean = mean(all_years_values, na.rm = TRUE), sd = sd(all_years_values, na.rm = TRUE))
    })
    names(climatic_mean_sd) <- climatic_pars # Name list elements by parameter


    # Function to normalize columns with mean and sd
    normalize <- function(train_data) {
        # Normalize non-climatic columns
        for (col in norm_cols) {
            # Standardize data
            standardized <- (train_data[[col]] - train_mean_sd[[paste0(col, "_mean")]]) /
                train_mean_sd[[paste0(col, "_sd")]]
            # Shift and scale to make positive
            standardized_min <- min(standardized, na.rm = TRUE)
            standardized_max <- max(standardized, na.rm = TRUE)
            train_data[[col]] <- (standardized - standardized_min) / (standardized_max - standardized_min)
        }

        # Normalize climatic columns using the global mean and sd for each parameter
        for (param in climatic_pars) {
            # Get column names for all years related to the parameter (e.g., srad_1985, srad_1986, etc.)
            param_cols <- grep(paste0("^", param, "_\\d{4}$"), names(train_data), value = TRUE)

            # Apply normalization for each column
            for (col in param_cols) {
                standardized <- (train_data[[col]] - climatic_mean_sd[[param]]$mean) /
                    climatic_mean_sd[[param]]$sd

                # Shift and scale to make positive
                standardized_min <- min(standardized, na.rm = TRUE)
                standardized_max <- max(standardized, na.rm = TRUE)
                train_data[[col]] <- (standardized - standardized_min) / (standardized_max - standardized_min)
            }
        }

        return(train_data)
    }

    train_data_norm <- normalize(train_data)
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


process_row <- function(cv_output, model, data_name, data_pars_name, biome_name, basic_pars_name) {

    # Define all possible parameters and identify missing ones
    all_possible_pars <- unique(unlist(c(non_data_pars, unique_colnames, categorical)))
    missing_cols <- setdiff(all_possible_pars, names(cv_output$pars))

    # Initialize row with provided and calculated values
    row <- data.frame(
        model = model,
        biome_name = biome_name,
        data = data_name,
        data_pars = data_pars_name,
        basic_pars = basic_pars_name,
        r2_mean = cv_output$r2_mean,
        r2_sd = cv_output$r2_sd,
        r2_final = cv_output$r2_final,
        stringsAsFactors = FALSE
    )

    # Add model parameters (coefficients or variable importance) from cv_output$pars
    pars_data <- as.data.frame(t(cv_output$pars))
    # Add missing columns with NA values
    for (col in missing_cols) {
        pars_data[[col]] <- NA
    }

    # Merge row with pars_data while ensuring consistent column order
    row <- cbind(row, pars_data[, all_possible_pars, drop = FALSE])

    # Define the desired order of columns and reorder
    desired_column_order <- c(
        "model", "biome_name", "data", "data_pars", "basic_pars",
        "r2_mean", "r2_sd", "r2_final", "age"
    )
    row <- row %>% select(all_of(desired_column_order), all_of(non_data_pars), everything())


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


find_combination_pars <- function(iter, data) {

    j <- iterations$data_par[iter]
    k <- iterations$biome[iter]
    l <- iterations$basic_par[iter]

    data_pars_iter <- data_pars[[k]][[j]]
    basic_pars_iter <- basic_pars[[l]]
    data <- normalize_independently(basic_pars_iter, data)$train_data

    # Initialize parameter vector with basic parameters and theta
    all_pars_iter <- c(setNames(
        rep(0, length(data_pars_iter)),
        c(data_pars_iter)
    ))

    all_pars_iter[["B0"]] <- mean(data[["biomass"]])
    all_pars_iter[["theta"]] <- 1
    basic_pars_iter <- c(basic_pars_iter, "theta")

    if ("age" %in% basic_pars_iter) {
        all_pars_iter["age"] <- 0
    }

    if ("k0" %in% basic_pars_iter) {
        all_pars_iter["k0"] <- -log(1 - mean(data[["biomass"]]) / mean(data[["nearest_mature_biomass"]]))
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
    best[["par"]] <- all_pars_iter[names(all_pars_iter) %in% basic_pars_iter]
    val <- 0

    base_row <- all_pars_iter
    base_row[names(all_pars_iter)] <- NA
    base_row <- c(likelihood = 0, base_row)

    should_continue <- TRUE
    # Iteratively add parameters and evaluate the model. Keep only AIC improvements.
    for (i in 1:length(data_pars_iter)) {
        if (!should_continue) break

        optim_remaining_pars <- foreach(j = remaining[-taken]) %dopar% {
        
        # iter_df <- data.frame()
        for (j in remaining[-taken]) {
            j = 20
            # check for categorical variables (to be included as a group)
            if (data_pars_iter[j] %in% categorical) {
                inipar <- c(best$par, all_pars_iter[grep(data_pars_iter[j], names(all_pars_iter))])
            } else {
                # as starting point, take the best values from last time
                inipar <- c(best$par, all_pars_iter[data_pars_iter[j]])
            }
            print(inipar)
            model <- run_optim(data, inipar, conditions)
            iter_row <- base_row
            iter_row[names(inipar)] <- model$par
            iter_row["likelihood"] <- model$value
            # return(iter_row)
        }

        iter_df <- as.data.frame(do.call(rbind, optim_remaining_pars))
        best_model <- which.min(iter_df$likelihood)
        # best_model_AIC <- 2 * iter_df$likelihood[best_model] + 2 * (i + val + 1)
        best_model_AIC <- iter_df$likelihood[best_model]

        if (best$AIC == 0 | best_model_AIC < best$AIC) {
            best$AIC <- best_model_AIC
            best$par <- iter_df[best_model, names(all_pars_iter)]
            best$par <- Filter(function(x) !is.na(x), best$par)
            taken <- which(sapply(data_pars_iter, function(x) any(grepl(x, names(best$par)))))
            print(paste0("iteration: ", iter, ", num parameters included: ", i, ", parameters taken: ", toString(data_pars_iter[taken])))
        } else {
            not_taken <- data_pars_iter[!data_pars_iter %in% names(best$par)]
            print(paste("No improvement. Exiting loop. Parameters not taken:", toString(not_taken)))
            should_continue <- FALSE
        }
    }

    return(best$par)
}
