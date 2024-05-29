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


  # - Imports the dataframe
  # - Removes unnecessary columns that will not be used in analysis
  # - Converts categorical data to dummy variables

import_data <- function(path, aggregate) {
  df <- read.csv(path)
  # Drop unnecessary columns
  df <- df[, -which(names(df) %in% c("system.index", ".geo", "biome"))]
  # create dummy variables for the categorical data with more than 2 types
  categorical <- c("ecoreg", "soil", "last_LU")
  df[categorical] <- lapply(df[categorical], as.factor)
  df <- createDummyFeatures(df, cols = categorical)
  if (aggregate == TRUE){
    # aggregate agbd by age
    df <- aggregate(agbd ~ age, df, median)
  }
  return(df)
}

# Intakes pars <- a vector with the initial parameters to be included
# data <- the dataframe with the predictors
# pars_chosen <- the list of parameters to be added into the shape term

# An example would look like:
#   growth_curve(c(B0 = 40, A = 80, theta = 5, age = 2), data_10[[1]], c('age'))

growth_curve <- function(pars, data, pars_chosen) {
  k <- data[[1]] * 0
  for (i in seq_along(pars_chosen)) {
    k <- k + pars[[pars_chosen[i]]] * data[[pars_chosen[i]]]
  }
  pars["B0"] + pars["A"] * (1 - exp(-k))^pars["theta"]
}


# Calculates Nonlinear Least Squares
# Intakes:
# fun <- the function to be used, either "nls" or "nll"
# pars <- a vector with the initial parameters to be included
# data <- the dataframe with the predictors
# pars_chosen <- the list of parameters to be added into the shape term
# conditions <- ranges of parameters to be restricted to

likelihood <- function(fun, pars, data, pars_chosen, conditions) {
  if (fun == "nll") {
    result <- -sum(dnorm(
      x = data$agbd - growth_curve(pars, data, pars_chosen), mean = 0,
      sd = pars["sd"], log = TRUE
    ), na.rm = TRUE)
    conditions <- c(conditions, 'pars["sd"] < 0')
  } else {
    result <- sum((growth_curve(pars, data, pars_chosen) - data$agbd)^2)
  }

  if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
    return(-Inf)
  } else if (is.na(result) || result == 0) {
    return(-Inf)
  } else {
    return(result)
  }
}

run_optimization <- function(fun, pars_basic, data, pars_chosen, conditions) {
  # Run optimization
  if (fun == "nll") {
    pars_basic <- c(pars_basic, setNames(0.1, "sd"))
  }
  o <- optim(c(pars_basic, setNames(0.1, pars_chosen[1])),
    likelihood,
    fun = fun,
    data = data,
    pars_chosen = pars_chosen[1],
    conditions = conditions
  )
  if (length(pars_chosen) > 1) {
    for (i in 2:length(pars_chosen)) {
      o <- optim(c(o$par, setNames(0.1, pars_chosen[i])),
        likelihood,
        fun = fun,
        data = data,
        pars_chosen = pars_chosen[1:i],
        conditions = conditions
      )
    }
  }
  print(paste("Residual sum of squares:", o$value))
  print(o$par)
  return(o)
}


################### Running model ###################

datafiles <- list(
  "data/unified_data_5_years.csv",
  "data/unified_data_10_years.csv",
  "data/unified_data_15_years.csv"
)

dataframes <- lapply(datafiles, import_data, aggregate = FALSE)
names_dataframes <- c("data_5", "data_10", "data_15")

# Define conditions
conditions <- list(
  'pars["theta"] > 10',
  'pars["theta"] < 0',
  'pars["B0"] < 0',
  'pars["B0"] > pars["A"]'
)

# intercept, asymptote, shape term, standard deviation
pars_basic <- c(B0 = 40, A = 80, theta = 5)

configurations <- list(
  c("age"), c("num_fires_before_regrowth"),
  c("age", "num_fires_before_regrowth"),
  c("age", "num_fires_before_regrowth", "all", "fallow", "indig", "protec"),
  setdiff(Reduce(intersect, lapply(dataframes, colnames)), c("b1", "agbd", "latitude", "longitude"))
)

names_configurations <- c("age", "fires", "age_fires", "all_cat", "all_non_LU", "all")

sum_squares_fit <- list()
pars_fit <- list()
# Run optimization
for (i in seq_along(configurations)) {
  for (j in seq_along(dataframes)) {
    print("----------------------------------------------------")
    print(names_dataframes[j])
    print(names_configurations[i])
    o_iter <- run_optimization(
      "nll", pars_basic, dataframes[[j]], configurations[[i]],
      if ("age" %in% configurations[[i]]) {
        c(conditions, list(
          'pars["age"] < 0',
          'pars["age"] > 5'
        ))
      } else {
        conditions
      }
    )
    sum_squares_fit[[paste(names_dataframes[j], names_configurations[i])]] <- o_iter$value
    pars_fit[[paste(names_dataframes[j], names_configurations[i])]] <- o_iter$par
  }
}

print(min(unlist(sum_squares)))
pars_fit[[15]]

pred <- growth_curve(pars_fit[[15]], dataframes[[3]], configurations[[5]])
plot(pred, dataframes[[3]]$agbd)
abline(0, 1)

# condit <- if ("age" %in% configurations[[i]]) {
#   c(conditions, list(
#     'pars["age"] < 0',
#     'pars["age"] > 5'
#   ))
# } else {
#   conditions
# }
