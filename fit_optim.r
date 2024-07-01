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

#------------------ SWITCHES ------------------#

run_all <- FALSE
run_one <- TRUE

#------------------ FUNCTIONS ------------------#

  # - Imports the dataframe
  # - Removes unnecessary columns that will not be used in analysis
  # - Converts categorical data to dummy variables

import_data <- function(path, aggregate) {
  df <- read.csv(path)
  # Drop unnecessary columns
  df <- df[, -which(names(df) %in% c("system.index", ".geo", "latitude", "longitude"))]
  # create dummy variables for the categorical data with more than 2 types
  categorical <- c("ecoreg", "soil", "last_LU")
  df[categorical] <- lapply(df[categorical], as.factor)
  df <- createDummyFeatures(df, cols = categorical)
  if (aggregate == TRUE) {
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
  pars["B0"] + (data["mat_gaus_ker"] - pars["B0"]) * (1 - exp(-k))^pars["theta"]
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
  # print(paste("Residual sum of squares:", o$value))
  # print(o$par)
  return(o)
}

make_asymptote <- function(df) {
  data <- read.csv("data/mature_biomass_climate_categ.csv")
  patterns <- c("si_", "prec_")
  means <- sapply(patterns, function(pat) rowMeans(data[, grep(pat, names(data))], na.rm = TRUE))
  colnames(means) <- c("mean_si", "mean_prec")
  data <- cbind(data, means)
  data <- data[, -grep("prec_|si_|biome|geo|system.index", names(data))]
  categorical <- c("ecoreg", "soil")
  data[categorical] <- lapply(data[categorical], as.factor)
  data <- createDummyFeatures(data, cols = categorical)
  data <- data %>%
    rename(agbd = b1, cwd = b1_1)

  # Identify numeric columns (excluding dummy variables)
  numeric_cols <- c("mean_si", "mean_prec", "cwd", "agbd")

  # Normalize numeric columns
  data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  })

  # Assuming `predicted_values` are your normalized predictions
  # And `min_x` and `max_x` are the stored minimum and maximum values for the column

  transform_back <- function(normalized_values, min_x, max_x) {
    original_values <- normalized_values * (max_x - min_x) + min_x
    return(original_values)
  }

  # Example usage
  # Let's say you have the min and max values for 'agbd' stored as `min_agbd` and `max_agbd`
  # And you have a vector of normalized predicted values for 'agbd' called `normalized_agbd`
  min_agbd <- min(data$agbd, na.rm = TRUE) # Assuming 'data' is your original dataset
  max_agbd <- max(data$agbd, na.rm = TRUE)

  original_agbd <- transform_back(normalized_agbd, min_agbd, max_agbd)

}

################### Running model ###################


datafiles <- list(
  "data/unified_data_5_years.csv",
  "data/unified_data_10_years.csv",
  "data/unified_data_15_years.csv"
)

dataframes <- lapply(datafiles, import_data, aggregate = FALSE)
names_dataframes <- c("data_5", "data_10", "data_15")
names(dataframes[[3]])



data <- dataframes[[3]]
# Get index of columns containing "prec" or "si"

# define the climatic parameters - the ones that change yearly
climatic_vars <- c("prec", "si")

# define the non-climatic parameters - the ones that are fixed throughout regrowth and
# that are used for fitting the model (excludes age and agbd)
non_climatic <- names(data)[!grepl("prec|si|agbd", names(data))]
clim_col_indices <- grep("prec", colnames(data))
tst <- max(clim_col_indices) + 1 - data["age"]

# list_of_lists <- list()
# for (i in seq_len(nrow(tst))) {
#   # Generate a sequence from the current tst element to max(clim_col_indices)
#   col_indices <- tst[i,]:max(clim_col_indices)
#   # Use the sequence to select the corresponding column names
#   col_names <- colnames(data)[col_indices]
#   # Add the column names to the list of lists
#   list_of_lists[[i]] <- col_names
# }

# pars <- c(prec = 0.1, si = 0.2)

# k <- data[[1]] * 0
# # Calculate the sum of climatic columns and non-climatic columns for each age
# for (age in 1:35) {
#   # Get the relevant years based on age
#   years <- seq(2019, 2019 - age + 1, by = -1)
#   clim_columns <- unlist(lapply(climatic_vars, function(pat) paste0(pat, "_", years)))
#   clim_columns
#   k <- lapply(climatic_vars, function(var) k + pars[var]) # data[[paste0(var, "_", years)]])
#   k
#   # Filter data for the current age
#   age_data <- data %>% filter(age == !!age)

#   for (i in c(1:age)) {
#     k <- k + pars[non_clim_var] * age_data[[non_clim_var]]
#   }

#   for (clim_var in climatic_var) {
#     k <- k + pars[clim_var] * data[[paste0(clim_var, "_", 2019 - age + 1)]]
#   }
# }

if (any(run_all, run_one)) {
  # Define conditions
  conditions <- list(
    'pars["theta"] > 10',
    'pars["theta"] < 0',
    'pars["B0"] < 0',
    'pars["B0"] > pars["A"]'
  )

  # intercept, asymptote, shape term, standard deviation
  pars_basic <- c(B0 = 40, A = 80, theta = 5)

  colnames_intersect <- Reduce(intersect, lapply(dataframes, colnames))
  colnames_filtered <- colnames_intersect[!grepl("b1|agbd|latitude|longitude|prec|si", colnames_intersect)]

  configurations <- list(
    c("age"),
    c("num_fires_before_regrowth"),
    c("age", "num_fires_before_regrowth"),
    c("age", "num_fires_before_regrowth", "all", "fallow", "indig", "protec"),
    colnames_filtered
  )
  names_configurations <- c("age", "fires", "age_fires", "all_cat", "all")
}

if (run_one) {
  pars_chosen <- configurations[[4]]
  pars <- c(
    pars_basic,
    setNames(rep(0.1, length(pars_chosen)), pars_chosen), setNames(0.1, "sd")
  )
  data = dataframes[[1]]
  result <- -sum(dnorm(
    x = data$agbd - growth_curve(pars, data, pars_chosen), mean = 0,
    sd = pars["sd"], log = TRUE
  ), na.rm = TRUE)

  run_optimization(
    "nll", pars_basic, data, pars_chosen,
    if ("age" %in% pars_chosen) {
      c(conditions, list(
        'pars["age"] < 0',
        'pars["age"] > 5'
      ))
    } else {
      conditions
    }
  )
}


if (run_all) {
  sum_squares_fit <- list()
  pars_fit <- list()
  # Run optimization
  for (i in seq_along(configurations)) {
    for (j in seq_along(dataframes)) {
      print("----------------------------------------------------")
      print(names_dataframes[j])
      print(names_configurations[i])
      o_iter <- run_optimization(
        "nls", pars_basic, dataframes[[j]], configurations[[i]],
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

  print(min(unlist(sum_squares_fit)))

  pred <- growth_curve(pars_fit[[15]], dataframes[[3]], configurations[[5]])
  plot(pred, dataframes[[3]]$agbd)
  abline(0, 1)

}

