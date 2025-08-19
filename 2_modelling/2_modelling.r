# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#
#                 Forest Regrowth Model Functions and Utilities
#
#                            Ana Avila - August 2024
#
#     This script defines the core functions used in the forest regrowth
#     modeling process (fit_3_run_model.r)
#
#     Functions included:
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------ Optimization for Forest Regrowth ------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
    if ("theta" %in% names(pars)) {
        conditions <- c(conditions, list('pars["theta"] < 0'))
    }

    return(optim(pars, calc_rss, data = train_data, conditions = conditions))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------- Sum of squared errors (SSE) ----------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   Calculates the sum of squared errors for the growth curve model,
#   incorporating parameter constraints.
#
# Arguments:
#   pars       : Named vector of parameter values to be evaluated.
#   data       : Dataframe containing predictor variables and forest attributes
#   conditions : List of parameter constraints expressed as character strings.
#
# Returns:
#   Numeric value representing the sum of squared errors.
#   Returns -Inf if constraints are violated or if the result is invalid.
#
# External Functions:
#   growth_curve()

calc_rss <- function(pars, data, conditions) {

    if ("lag" %in% names(pars)) {
        growth_curves <- growth_curve(pars, data, pars[["lag"]])
        residuals <- growth_curves - data$biomass
        result <- sum(residuals^2)
    } else {
        result <- sum((growth_curve(pars, data) - data$biomass)^2)
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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------- Chapman-Richards Growth Curve -------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#
# Function Description:
#   Calculates the Chapman-Richards growth curve based on provided parameters and data.
#   Incorporates intercept terms and forest age.

# get the mean value per column of all columns in data with the name pdsi
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



growth_curve <- function(pars, data, lag = 0) {

    fit_data_pars <- setdiff(names(pars), c(non_data_pars, "age", "theta"))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k
    k <- rep(pars[["k0"]], nrow(data))

    age <- data[["age"]]

    if ("lag" %in% names(pars)) {
        pars[["B0"]] <- 0
        age <- age + lag
    }

    if (length(fit_data_pars) > 0) {
        k <- (k + rowSums(sapply(fit_data_pars, function(par) {
            pars[[par]] * data[[par]]
        }, simplify = TRUE))) * (age)
    } else {
        k <- k * age
    }

    # Constrains k to avoid negative values
    k[which(k < 1e-10)] <- 1e-10
    k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)

    if ("theta" %in% names(pars)) {
        theta <- pars[["theta"]]
        pars[["B0"]] <- 0
    } else {
        theta <- 1
    }

    return(pars[["B0"]] + (data[["asymptote"]] - pars[["B0"]]) * (1 - exp(-k))^theta)
}

