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

    if ("B0" %in% names(pars)) {
        conditions <- c(conditions, list('pars["B0"] < 0'))
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
    # Define parameters that are not expected to change yearly (not prec or si)
    non_clim_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k
    if ("lag" %in% names(pars)) {
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
            # as.matrix(t()) is used to ensure that rowSums would work in cases with a single row
            k[indices] <- k[indices] + rowSums(as.matrix(t(sapply(clim_columns, function(col) pars[[clim_par]] * data[[col]][indices]))))
        }
    }

    # Constrains k to avoid negative values
    k[which(k < 1e-10)] <- 1e-10
    k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)

    return(pars[["B0"]] + (data[["nearest_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
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
        # Calculate log-normal scaled base using re_base and parameters
        growth_curves <- growth_curve(pars, data, exp(pars["lag"]))
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
#   data             : The dataframe to be used for analysis
# Returns:
#   data             : A dataframe with normalized numerical values

normalize_independently <- function(train_data, test_data = NULL) {

    # Select numeric columns for normalization, excluding specified ones
    exclusion_list <- c(unlist(categorical), "biomass", "nearest_biomass")
    # if k is multiplied by the age column, don't normalize age
    exclusion_list <- c(exclusion_list, if (!"age" %in% names(pars)) "age")

    norm_cols <- c(names(train_data)[!grepl(paste0(exclusion_list, collapse = "|"), names(train_data))])

    # Compute mean and standard deviation for normalization based on training data
    train_mean_sd <- train_data %>%
        summarise(across(all_of(norm_cols), list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE))))

    if (any(climatic_pars %in% names(data))){
        # Compute global mean and standard deviation for each climatic parameter across all years
        climatic_mean_sd <- lapply(climatic_pars, function(param) {
            all_years_values <- train_data %>%
                select(matches(paste0("^", param, "_\\d{4}$"))) %>%
                unlist(use.names = FALSE) # Flatten to a single vector
            list(mean = mean(all_years_values, na.rm = TRUE), sd = sd(all_years_values, na.rm = TRUE))
        })
        names(climatic_mean_sd) <- climatic_pars # Name list elements by parameter
    }

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
        
        if (any(climatic_pars %in% names(data))) {
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

    all_pars[["theta"]] <- 1
    all_pars["k0"] <- -log(1 - mean(data[["biomass"]]) / mean(data[["nearest_biomass"]]))

    if ("lag" %in% basic_pars) {
        all_pars["lag"] <- 0
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
            # check for categorical variables (to be included as a group)
            if (data_pars[j] %in% categorical) {
                inipar <- c(best$par, all_pars[grep(data_pars[j], names(all_pars))])
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
