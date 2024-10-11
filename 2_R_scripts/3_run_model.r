library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)

# Source external R scripts for data import and function definitions
source("./2_R_scripts/1_import_data.r")
source("./2_R_scripts/2_functions.r")

re_base <- dnorm(1000)

# Import data
# data <- read_csv("./0_data/eu.csv", show_col_types = FALSE)
# data <- read_csv("./0_data/non_aggregated.csv", show_col_types = FALSE)
data <- import_data("./0_data/non_aggregated.csv", convert_to_dummy = TRUE)
data <- data %>% filter(biome == 4)

# df <- subset(data, select = -c(agbd, biome))
# par_columns <- colnames(df)
par_columns <- c("sur_cover")

basic_pars <- c("k0", "m_base", "sd_base", "theta", "sd")

pars <- setNames(rep(0.0001, length(par_columns)), par_columns)
pars[["theta"]] <- 1
pars[["k0"]] <- 0
pars[["m_base"]] <- 0
pars[["sd_base"]] <- 1
pars[["sd"]] <- 1

conditions <- list(
    function(pars) pars["theta"] > 10,
    function(pars) pars["theta"] < 0,
    function(pars) pars["sd"] < 0,
    function(pars) pars["sd_base"] < 0,
    function(pars) pars["m_base"] < 0
)

growth_curve_lag <- function(pars, data, lag) {
    k <- rep(pars["k0"], nrow(data))

    # Get the names of data_pars to iterate over
    data_pars_names <- names(pars)[!names(pars) %in% basic_pars]

    if (length(data_pars_names) > 0) {
        k <- k + rowSums(sapply(data_pars_names, function(par) {
            pars[[par]] * data[[par]]
        }, simplify = TRUE)) * (data[["age"]] + lag)
    }

    return(data[["nearest_mature_biomass"]] * (1 - exp(-k))^pars[["theta"]])
}

likelihood_lag <- function(pars, data, conditions) {

    # Calculate log-normal scaled base using re_base and parameters
    scaled_base <- exp((re_base + pars["m_base"]) * pars["sd_base"])

    m_results <- matrix(0, nrow = nrow(data), ncol = length(scaled_base))

    growth_curves <- sapply(scaled_base, function(lag) growth_curve_lag(pars, data, lag))
    differences <- sweep(growth_curves, 1, data$agbd, "-")
    m_results <- dnorm(differences, sd = pars["sd"])

    mean_results <- rowMeans(m_results)
    mean_results[mean_results == 0] <- 1e-10

    # Calculate log-likelihood
    result <- sum(-log(mean_results))

    # Check parameter constraints
    if (any(sapply(conditions, function(cond) cond(pars)))) {
        return(-Inf)
    } else if (is.na(result) || result == 0) {
        return(-Inf)
    } else {
        return(result)
    }
}

# Run the optimization using the optim function
model <- optim(par = pars, fn = likelihood_lag, data = data, conditions = conditions)

calc_rsq <- function(data, pred) {
    obs_pred <- lm(data$agbd ~ pred)
    residuals <- summary(obs_pred)$residuals
    sum_res_squared <- sum(residuals^2)
    total_sum_squares <- sum((data$agbd - mean(data$agbd))^2)
    rsq <- 1 - (sum_res_squared / total_sum_squares)

    return(rsq)
}

calc_rsq(data, unlist(growth_curve(model$par, data, exp((re_base + model$par["m_base"]) * model$par["sd_base"]))))


growth_curve_B0_theta <- function(pars, data) {
    k <- rep(pars["k0"], nrow(data))

    # Get the names of data_pars to iterate over
    data_pars_names <- names(pars)[!names(pars) %in% basic_pars]

    if (length(data_pars_names) > 0) {
        k <- k + rowSums(sapply(data_pars_names, function(par) {
            pars[[par]] * data[[par]]
        }, simplify = TRUE)) * (data[["age"]])
    }

    return(pars[["B0"]] + (data[["nearest_mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
}

likelihood_B0_theta <- function(pars, data, conditions) {
    # Calculate log-normal scaled base using re_base and parameters
    scaled_base <- exp((re_base + pars["m_base"]) * pars["sd_base"])

    m_results <- matrix(0, nrow = nrow(data), ncol = length(scaled_base))

    growth_curves <- sapply(scaled_base, function(lag) growth_curve_lag(pars, data, lag))
    differences <- sweep(growth_curves, 1, data$agbd, "-")
    m_results <- dnorm(differences, sd = pars["sd"])

    mean_results <- rowMeans(m_results)
    mean_results[mean_results == 0] <- 1e-10

    # Calculate log-likelihood
    result <- sum(-log(mean_results))

    # Check parameter constraints
    if (any(sapply(conditions, function(cond) cond(pars)))) {
        return(-Inf)
    } else if (is.na(result) || result == 0) {
        return(-Inf)
    } else {
        return(result)
    }
}

par_columns <- c("sur_cover")

basic_pars <- c("k0", "B0", "theta", "sd")

pars <- setNames(rep(0.0001, length(par_columns)), par_columns)
pars[["theta"]] <- 1
pars[["k0"]] <- 0
pars[["B0"]] <- mean(data$agbd)