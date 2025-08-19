# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#        Forest Regrowth Model Functions and Utilities
#
#                   Ana Avila - August 2024
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------ Optimization for Forest Regrowth ------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Optim wrapper function.
#' Applies constraints and optimizes the Chapman-Richards
#' growth model to forest regrowth data.
#'
#' @param train_data Data frame. Training dataset containing forest attributes
#'   such as age, biomass, and predictors.
#' @param pars Named numeric vector. Initial parameter values 
#' (e.g., B0, k0, theta, lag).
#' @param conditions List of character strings. Parameter constraints expressed
#'   as logical conditions (evaluated during optimization).
#'
#' @return A list. The output object from \code{optim()}, including
#' optimized parameters, objective value, and convergence status.


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
# --------------- Sum of squared errors (SSE) --------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Computes the residual sum of squares between observed biomass values
#' and predictions from the Chapman-Richards growth curve.
#'
#' @param pars Named numeric vector. Candidate parameter values.
#' @param data Data frame. Must include \code{age}, \code{biomass}, 
#' \code{asymptote}, and predictors as needed.
#' @param conditions List of character strings. Constraints on parameters 
#' expressed as logical conditions.
#'
#' @return Numeric. The residual sum of squares value, or \code{-Inf}
#' if constraints are violated or if residuals are invalid.
#'
#' @details
#' - Constraint violations or invalid objective values automatically 
#' return \code{-Inf}, so that they are rejected during optimization.


calc_rss <- function(pars, data, conditions) {

    if ("lag" %in% names(pars)) {
        growth_curves <- growth_curve(pars, data, pars[["lag"]])
        residuals <- growth_curves - data$biomass
        result <- sum(residuals^2)
    } else {
        result <- sum((growth_curve(pars, data) - data$biomass)^2)
    }

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

#' Estimates forest biomass according to the Chapman-Richards equation,
#' with optional predictor effects and lag adjustment.
#'
#' @param pars Named numeric vector. Growth model parameters:
#'   - \code{B0} : Initial biomass at age zero (baseline).
#'   - \code{k0} : Baseline growth rate constant.
#'   - \code{theta} : Curve shape parameter.
#'   - \code{lag} : Optional regrowth lag (age offset).
#'   - Additional coefficients for predictor covariates.
#' @param data Data frame. Must include \code{age}, \code{biomass}, \code{asymptote}, and predictors if present.
#' @param lag Numeric (default = 0). Optional time adjustment for regrowth start.
#'
#' @return Numeric vector. Predicted biomass values for each row in \code{data}.
#'
#' @details
#' - Growth rate (\code{k}) is computed from \code{k0}, covariates, and forest age.
#' - Growth rates are constrained for numerical stability:
#'   values below \code{1e-10} are floored, and values above 7 are capped.
#' - If \code{theta} is specified, \code{B0} is internally fixed at 0 to avoid redundancy.
#'

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
    # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)
    k[which(k > 7)] <- 7

    if ("theta" %in% names(pars)) {
        theta <- pars[["theta"]]
        pars[["B0"]] <- 0
    } else {
        theta <- 2
    }

    return(pars[["B0"]] + (data[["asymptote"]] - pars[["B0"]]) * (1 - exp(-k))^theta)
}

