####################################################################
############### Process data and define functions ##################
####################################################################
# Ana Avila - May 2024
# ~~~~~~~~~~~~~~~~~~~~
# Intakes:
# Outputs:
####################################################################

# get all R squareds and save them properly
# include nearest neighbor distances for wherever the kernel gives zero
# how many of the secondary forests distance histogram
# fit the kernel with sigma in R for the gaussian kernel

library(ggplot2)
library(terra)
library(tidyverse) # for stringr operations
library(mlr) # for createDummyFeatures

#------------------ SWITCHES ------------------#

run_all <- FALSE
run_one <- TRUE
include_k_par <- TRUE

#------------------ FUNCTIONS ------------------#

# - Imports the dataframe
# - Removes unnecessary columns that will not be used in analysis
# - Converts categorical data to dummy variables

import_data <- function(path, aggregate, normalize) {
  df <- read.csv(path)
  # Drop unnecessary columns
  df <- df[, -which(names(df) %in% c("system.index", ".geo", "latitude", "longitude", "biome"))]
  # create dummy variables for the categorical data with more than 2 types
  categorical <- c("ecoreg", "soil", "last_LU")
  df[categorical] <- lapply(df[categorical], as.factor)
  df <- createDummyFeatures(df, cols = categorical)
  if (aggregate == TRUE) {
    # aggregate agbd by age
    df <- aggregate(agbd ~ age, df, median)
  }

  df_climatic_hist <- data.frame()
  # keep only the climatic variables of the years in which there was regrowth
  for (age in 1:max(df$age)) {
    age_data <- df[df$age == age, ]

    # Get the relevant years based on age
    years <- seq(2019, 2019 - age + 1, by = -1)
    clim_columns <- unlist(lapply(climatic_vars, function(pat) paste0(pat, "_", years)))

    # Identify all columns including "prec_" or "si_"
    # Assuming climatic_vars is defined as
    # climatic_vars <- c("prec", "si")

    # Create a regex pattern from climatic_vars to match column names
    pattern <- paste(climatic_vars, collapse = "|")
    pattern <- paste0("(", pattern, ")_")

    # Use the pattern in grep to find all climatic columns
    all_clim_columns <- grep(pattern, colnames(df), value = TRUE)

    # Exclude columns that are in clim_columns
    clim_columns_not_included <- all_clim_columns[!all_clim_columns %in% clim_columns]

    age_data[clim_columns_not_included] <- lapply(age_data[clim_columns_not_included], function(x) 0)

    df_climatic_hist <- rbind(df_climatic_hist, age_data)
  }
  # data <- data[, -which(names(data) == "lulc_sum_62")]
  if (normalize) {
    # Normalize continuous variables
    continuous_cols <- names(df_climatic_hist)[
      !grepl("soil|ecoreg|last_LU|protec|indig|agbd|mat_gaus_ker", names(df_climatic_hist))
    ]
    df_climatic_hist[continuous_cols] <- lapply(df_climatic_hist[continuous_cols], function(x) {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    })
    df_climatic_hist <- df_climatic_hist[, colSums(is.na(df_climatic_hist)) < nrow(df_climatic_hist)]
  }

  return(df_climatic_hist)
}

# Intakes pars <- a vector with the initial parameters to be included
# data <- the dataframe with the predictors
# pars_chosen <- the list of parameters to be added into the shape term

# An example would look like:
#   growth_curve(c(B0 = 40, A = 80, theta = 5, age = 2), data_10[[1]], c('age'))

# define the climatic parameters - the ones that change yearly
growth_curve <- function(pars, data, pars_chosen) {
  k <- data[[1]] * 0 
  if (include_k_par) {
    k <- k + pars[["k"]]
  }
  non_clim_vars <- pars_chosen[!pars_chosen %in% climatic_vars]

  if (any(climatic_vars %in% pars_chosen)) {
      k <- k + rowSums(sapply(non_clim_vars, function(var) pars[[var]] * data[[var]] * data[["age"]]))

    for (clim_var in climatic_vars) {
      years <- seq(2019, 1985, by = -1)
      clim_columns <- unlist(paste0(clim_var, "_", years))

      k <- k + rowSums(sapply(clim_columns, function(col) pars[[clim_var]] * data[[col]]))
    }
  } else {
    k <- k + rowSums(sapply(non_clim_vars, function(var) pars[[var]] * data[[var]]))
  }

  pars[["B0"]] + (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]]
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
  return(o)
}

################### Running model ###################

datafiles <- list(
  "data/15y_LULC.csv"
)

climatic_vars <- c("prec", "si")
dataframes <- lapply(datafiles, import_data, aggregate = FALSE, normalize = FALSE)
# names_dataframes <- c("data_5", "data_15")
data <- dataframes[[1]]
data <- data[!is.na(data$mature_biomass), ]

if (any(run_all, run_one)) {
  # Define conditions
  conditions <- list(
    'pars["theta"] > 10',
    'pars["theta"] < 0',
    'pars["B0"] < 0'
  )

  # intercept, shape term, growth rate
  pars_basic <- c(B0 = 40, theta = 5)
  if (include_k_par){
    pars_basic <- c(pars_basic, k = 0.01)
  }

  colnames_intersect <- Reduce(intersect, lapply(dataframes, colnames))
  colnames_filtered <- colnames_intersect[!grepl("agbd|latitude|longitude|prec|si|mat_gaus_ker", colnames_intersect)]

  configurations <- list(
    c("age"),
    c("num_fires_before_regrowth"),
    c("age", "num_fires_before_regrowth"),
    c("age", "num_fires_before_regrowth", "sur_cover", "fallow", "indig", "protec", "cwd"),
    colnames_filtered,
    c(colnames_filtered, climatic_vars)
  )
  names_configurations <- c("age", "fires", "age_fires", "all_cat", "all_non_hist", "all_hist")
}


if (run_one) {
  pars_chosen <- configurations[[3]]
  pars <- c(
    pars_basic,
    setNames(rep(0.1, length(pars_chosen)), pars_chosen)#, setNames(0.1, "sd")
  )
  # data = dataframes[[4]]
  # result <- -sum(dnorm(
  #   x = data$agbd - growth_curve(pars, data, pars_chosen), mean = 0,
  #   sd = pars["sd"], log = TRUE
  # ), na.rm = TRUE)

  val <- run_optimization(
    "nls", pars_basic, data, pars_chosen,
    if ("age" %in% pars_chosen) {
      c(conditions, list(
        'pars["age"] < 0',
        'pars["age"] > 5'
      ))
    } else {
      conditions
    }
  )

  pred <- growth_curve(val$par, data, pars_chosen)
  rsq_fun <- function(x, y) cor(x, y)^2
  rsq_fun(data$agbd, pred)
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
      pred <- growth_curve(o_iter$par, dataframes[[j]], configurations[[i]])
      rsq <- rsq_fun(dataframes[[j]]$agbd, pred)
      print(rsq)
      sum_squares_fit[[paste(names_dataframes[j], names_configurations[i])]] <- rsq
      pars_fit[[paste(names_dataframes[j], names_configurations[i])]] <- o_iter$par
    }
  }

  print(min(unlist(sum_squares_fit)))

  # pred <- growth_curve(pars_fit[[15]], dataframes[[3]], configurations[[5]])
  # plot(pred, dataframes[[3]]$agbd)
  # abline(0, 1)

}