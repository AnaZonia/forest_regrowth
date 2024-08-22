# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------ Growth Curve -------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#     This function calculates the Chapman-Richards growth curve based on the provided
#     parameters and data. It can incorporate:
#     - Yearly-changing climatic parameters (prec and si)
#     - The intercept term defined either as B0 or B0_exp (out or in of the exponential)
#     - A growth rate intercept term k0
#     - Ages of secondary forests as an explicit parameter (part of "pars")
#       or as an implicit parameter (multiplying all non-yearly predictors by age)
#
# Intakes:
#   pars <- a vector with the named parameters to be included
#   data <- the chosen training dataframe


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

    k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)

    if (fit_logistic) {
        # return(data[["mature_biomass"]] * (1 / (1 + exp(k))))
        return(pars[["B0"]] * data[["mature_biomass"]] * exp(k)) / ((data[["mature_biomass"]] - B0) + B0 * exp(k))
    } else {
        if ("B0" %in% names(pars)) {
            return(pars[["B0"]] + (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
        } else {
            return(data[["mature_biomass"]] * (1 - exp(-(pars[["B0_exp"]] + k)))^pars[["theta"]])
        }
    }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------------------------- Likelihood -------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Intakes:
#   fun <- the function to be used, either "nls" (Nonlinear Least Squares)
#          or "nll" (Negative Log Likelihood)
#   pars <- a vector with the named parameters to be included
#   data <- the chosen training dataframe
#   conditions <- ranges of parameters to be restricted to

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
# ---------------------------------- Optimization ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#     This function prepares the optimization process for the forest regrowth model.
#     It runs the optimization process iteratively, incorporating one parameter of the
#     chosen set at a time.
#
#     It also determines:
#     - A new "sd" parameter if the function is set to "nll" (Negative Log Likelihood)
#     - The theta ("shape") parameter of the Chapman-Richards growth curve
#     - The values of the initial parameters (0.1 for all except for theta and B0)
#
# Intakes:
#   fun <- the function to be used, either "nls" (Nonlinear Least Squares)
#          or "nll" (Negative Log Likelihood)
#   basic_pars_iter <- a vector with the named parameters that do not correspond to data,
#                 but are to be included as part of the functional form
#   data_pars_iter <- a vector with the chosen parameter set to be included
#                  (corresponding to columns in dataframe)
#   train_data <- the chosen training dataframe
#   test_data <- the chosen testing dataframe (for rsq calculation)
#   conditions <- ranges of parameters to be restricted to

run_optim <- function(train_data, test_data, pars, conditions) {
    if ("age" %in% names(pars)) {
        conditions <- c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
    }
    if ("B0" %in% names(pars)) {
        conditions <- c(conditions, list('pars["B0"] < 0'))
    }

    model <- optim(pars, likelihood, data = train_data, conditions = conditions)

    filtered_test_data <- filtered_data(train_data, test_data)
    pred <- growth_curve(model$par, filtered_test_data)
    rsq <- calc_rsq(filtered_test_data, pred)
    print(paste("R-squared:", rsq))

    return(list(
        model_par = model$par,
        rsq = rsq
    ))
}



run_lm <- function(train_data, test_data, pars) {
    lm_formula <- as.formula(paste("agbd ~", paste(pars, collapse = " + ")))

    model <- lm(lm_formula, data = train_data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Check for rank deficiency, and remove variables if necessary
    # (usually ones with too few unique values)
    aliased_vars <- summary(model)$aliased

    if (any(aliased_vars)) {
        problematic_vars <- names(aliased_vars)[aliased_vars]
        print(paste("Removing rank-deficient variables:", paste(problematic_vars, collapse = ", ")))

        # Remove problematic variables from the formula
        pars <- pars[!pars %in% problematic_vars]
        lm_formula <- as.formula(paste("agbd ~", paste(pars, collapse = " + ")))
        model <- lm(lm_formula, data = train_data)
    }

    filtered_test_data <- filtered_data(train_data, test_data)
    pred <- predict(model, newdata = filtered_test_data)
    rsq <- calc_rsq(filtered_test_data, pred)
    print(paste("R-squared:", rsq))

    return(list(
        model_par = t(summary(model)$coefficients[-1, 1, drop = FALSE]), # -1 to remove (Intercept),
        rsq = rsq
    ))
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----------------- Ensuring testing data is within testing data range -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

filtered_data <- function(train_data, test_data) {
    # Calculate min and max for each column in train_data
    train_min <- sapply(train_data, min)
    train_max <- sapply(train_data, max)

    # Function to check if a row is within the range
    is_within_range <- function(row) {
        all(row >= train_min & row <= train_max)
    }

    # Apply the function to each row of test_data
    filtered_test_data <- test_data[apply(test_data, 1, is_within_range), ]

    return(filtered_test_data)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----------------- Calculating R squared -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

calc_rsq <- function(data, pred) {
    obs_pred <- lm(data$agbd ~ pred)
    residuals <- summary(obs_pred)$residuals
    sum_res_squared <- sum(residuals^2)
    total_sum_squares <- sum((data$agbd - mean(data$agbd))^2)
    rsq <- 1 - (sum_res_squared / total_sum_squares)

    return(rsq)
}
