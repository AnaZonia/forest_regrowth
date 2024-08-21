library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)

# Source external R scripts for data import and function definitions
source("fit_1_import_data.r")
source("fit_2_functions.r")

set.seed(1)

# Define land-use history intervals to import four dataframes
intervals <- list("5yr", "10yr", "15yr", "all")
datafiles_amaz <- paste0("./data/amaz_", intervals, ".csv")
dataframes_amaz <- lapply(datafiles_amaz, import_climatic_data, normalize = TRUE)
datafiles_countrywide <- paste0("./data/countrywide_", intervals, ".csv")
dataframes_countrywide <- lapply(datafiles_countrywide, import_climatic_data, normalize = TRUE)

run_growth_model <- function(data, initial_pars) {
    conditions <- list(
        'pars["theta"] > 10'
        ,'pars["theta"] < 0'
        ,'pars["B0"] < 0'
        ,'pars["age"] < 0'
        ,'pars["age"] > 5'
    )

    k <- c(pars[["age"]] * data[["age"]] + pars[["cwd"]] * data[["cwd"]])

    growth_curve <- function(pars, data) {
        pars[["B0"]] + (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]]
    }

    likelihood <- function(pars, data, conditions) {
        result <- sum((growth_curve(pars, data) - data$agbd)^2)

        if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
            return(Inf)
        } else if (is.na(result) || result == 0) {
            return(Inf)
        } else {
            return(result)
        }
    }

    o <- optim(initial_pars, likelihood, data = data, conditions = conditions)

    pred <- growth_curve(o$par, data)
    rsq <- cor(data$agbd, pred)^2

    return(list(
        optimized_parameters = o$par,
        r_squared = rsq,
        predictions = pred,
        optimization_result = o
    ))
}


for (i in seq_along(dataframes_amaz)) {
  print("----------------------------------------------------")
  print(i)

  # Without asymptote (using mature_biomass)
  result1 <- run_growth_model(dataframes_amaz[[i]], c(B0 = 40, theta = 5, age = 0, cwd = 0))
  print(paste("R-squared amaz:", result1$r_squared))

  # # With asymptote
  result2 <- run_growth_model(dataframes_countrywide[[i]], c(B0 = 40, theta = 5, age = 0, cwd = 0))
  print(paste("R-squared countrywide_amaz_subset:", result2$r_squared))
}
