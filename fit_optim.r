####################################################################
############### Process data and define functions ##################
####################################################################
# Ana Avila - May 2024
# ~~~~~~~~~~~~~~~~~~~~
# Intakes:
# Outputs:
####################################################################

library(ggplot2)
library(terra)
library(tidyverse) # for stringr operations
library(mlr) # for createDummyFeatures

"
  - Imports the dataframe
  - Removes unnecessary columns that will not be used in analysis
  - Converts categorical data to dummy variables
"
import_data <- function(path) {
  df <- read.csv(path)

  # Drop unnecessary columns
  df <- df[, -which(names(df) %in% c("system.index", ".geo", "biome"))]
  # create dummy variables for the categorical data with more than 2 types
  categorical <- c("ecoreg", "soil", "last_LU")
  df[categorical] <- lapply(df[categorical], as.factor)
  df <- createDummyFeatures(df, cols = categorical)
  # aggregate agbd by age
  df_aggregated <- aggregate(agbd ~ age, df, median)

  return(list(data = df, data_aggregated = df_aggregated))
}

growth_curve <- function(pars, data, pars_chosen) {
  k <- 0
  for (i in seq_along(pars_chosen)) {
    k <- k + pars[[pars_chosen[i]]] * data[[pars_chosen[i]]]
  }
  pars["B0"] + pars["A"] * (1 - exp(-k))^pars["theta"]
}

write_init_parameters <- function(pars_basic, pars_chosen) {
  pars_init_values <- rep(0.1, length(pars_chosen))
  new_elements <- setNames(pars_init_values, pars_chosen)
  return(c(pars_basic, new_elements))
}

# Nonlinear Least Squares
nls <- function(pars, data, growth_curve, conditions, pars_chosen) {
  for (condition in conditions) {
    if (!eval(condition)) {
      return(-Inf)
    }
  }
  result <- sum((growth_curve(pars, data, pars_chosen) - data$agbd)^2)
  if (is.na(result) || result == 0) {
    return(-Inf)
  } else {
    return(result)
  }
}

fit_optim <- function(pars, fn, data, growth_curve, conditions) {
  o <- optim(pars, fn, data = data, growth_curve = growth_curve, conditions = conditions)
  print(paste("Residual sum of squares:", o$value))
  print(o$par)
  return(o)
}

################### Running model ###################

data_5 <- import_data("data/unified_data_5_years.csv")
data_10 <- import_data("data/unified_data_10_years.csv")
data_15 <- import_data("data/unified_data_15_years.csv")

conditions <- list(
  expression(pars["age"] < 0),
  expression(pars["age"] > 5),
  expression(pars["theta"] > 10),
  expression(pars["theta"] < 0),
  expression(pars["B0"] < 0)
)

# intercept, asymptote, shape term, standard deviation
pars_basic <- c(B0 = 40, A = 80, theta = 5)
age_cwd <- c("age", "b1")
age <- c("age")

# Run optimization
pred_5_age <- fit_optim(
  pars = write_init_parameters(pars_basic, age),
  nls,
  data = data_5[[1]],
  growth_curve = growth_curve,
  conditions = conditions,
  pars_chosen = age
)

growth_curve(write_init_parameters(pars_basic, age), data_5[[1]], age)
nls(write_init_parameters(pars_basic, age), data_5[[1]], growth_curve, conditions, age)
result <- sum((growth_curve(write_init_parameters(pars_basic, age), data_5[[1]], age) - data_5[[1]]$agbd)^2)

# pred <- fit_optim(write_init_parameters(pars_basic, 2), nls, data = data_5, growth_curve = function(pars, data) {
#   growth_curve(pars, data, c("age", "num_fires_before_regrowth"))
# })

# pred_simple_agg <- fit_optim(pars_simple, nls,
#   data = data_aggregated,
#   growth_curve = growth_curve_1_pred
# )
# pred_simple <- fit_optim(pars_simple, nls,
#   data = data, growth_curve = growth_curve_1_pred
# )




# pred_numfires <- fit_optim(pars_env_par, nls, data = data,
# growth_curve = function(pars, data) {
#   growth_curve_1_pred(pars, data, "num_fires_before_regrowth")
# })
# pred_sur_for_cov <- fit_optim(pars_env_par, nls, data = data, growth_curve = function(pars, data) {
#   growth_curve_1_pred(pars, data, "all")
# })

# pred_simple$value - pred_numfires$value
# pred_simple$value - pred_sur_for_cov$value

# pred <- growth_curve(o$par, data_aggregated)

# plot(data$agbd, pred)
# abline(c(0, 1))