

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
    if ("theta" %in% names(pars)) {
        conditions <- c(conditions, list('pars["theta"] < 0'))
    }

    return(optim(pars, calc_rss, data = train_data, conditions = conditions))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------ Likelihood Function for Optimization -------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------- Chapman-Richards Growth Curve ----------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function Description:
#   Calculates the Chapman-Richards growth curve based on provided parameters and data.
#   Incorporates yearly-changing climatic parameters, intercept terms, and forest age.

# get the mean value per column of all columns in data with the name pdsi
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


growth_curve <- function(pars, data, lag = 0) {

    # I am not checking climatic variables correctly - need to add monthly or yearly separately.
    # Define parameters that are not expected to change yearly (not prec or si)
    non_yearly_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars, "age"))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the growth rate k
    k <- rep(pars[["k0"]], nrow(data))

    age <- data[["age"]]

    if ("lag" %in% names(pars)) {
        pars[["B0"]] <- 0
        age <- age + lag
    }

    if (length(non_yearly_pars) > 0) {
        k <- (k + rowSums(sapply(non_yearly_pars, function(par) {pars[[par]] * data[[par]]}, simplify = TRUE))) * (age)
    } else {
        k <- k * age
    }


    # Add yearly-changing climatic parameters to the growth rate k (if included in the parameter set)
    for (clim_par in intersect(climatic_pars, names(pars))) {
        last_year <- 1985
        year_seq <- seq(last_year, 2019)
        clim_columns <- paste0(clim_par, "_", year_seq)
        # print(clim_columns)

        k[indices] <- k[indices] + pars[[clim_par]] * rowMeans(sapply(clim_columns, function(col) data[[col]][indices]))

        # for (yr in 1:max(data[["age"]])) {

        #     indices <- which(data[["age"]] == yr)
        #     # Generate a sequence of years for the current age group
        #     # Starting from 2019 and going back 'yrs' number of years (stopping in 1958 as the earliest year)
        #     last_year <- max(2019 - yr - round(lag) + 1, 1959)
        #     lag_rounded <- max(round(lag), 0)
        #     last_year <- max(2019 - yr - lag_rounded + 1, 1959)
        #     year_seq <- seq(last_year, 2019)
        #     clim_columns <- paste0(clim_par, "_", year_seq)

        #     # Add the contribution individually or as an average

        #     k[indices] <- k[indices] + pars[[clim_par]] * rowMeans(sapply(clim_columns, function(col) data[[col]][indices]))

        #     # k[indices] <- k[indices] + pars[[clim_par]] * rowSums(sapply(clim_columns, function(col) data[[col]][indices]))
        # }
    }

    # Constrains k to avoid negative values
    k[which(k < 1e-10)] <- 1e-10
    k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)

    if ("theta" %in% names(pars)) {
        theta <- pars[["theta"]]
    } else {
        theta <- 1.1
    }

    return(pars[["B0"]] + (data[["asymptote"]] - pars[["B0"]]) * (1 - exp(-k))^theta)
}

