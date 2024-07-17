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
library(foreach)
library(doParallel)

source("1_import_data.r")

#------------------ SWITCHES ------------------#

run_all <- TRUE
run_one <- TRUE

#------------------ FUNCTIONS ------------------#


# Intakes pars <- a vector with the initial parameters to be included
# data <- the dataframe with the predictors
# pars_chosen <- the list of parameters to be added into the shape term

# define the climatic parameters - the ones that change yearly
growth_curve <- function(pars, data, pars_chosen) {
  if ("k0" %in% pars_chosen) {
    k <- rep(pars[["k0"]], nrow(data))
    non_clim_vars <- setdiff(pars_chosen, c("k0", climatic_vars))
  } else {
    k <- rep(0, nrow(data))
    non_clim_vars <- setdiff(pars_chosen, climatic_vars)
  }

  implicit_age <- 1
  if ("implicit_age" %in% pars_chosen) {
    implicit_age <- data[["age"]]
    non_clim_vars <- setdiff(non_clim_vars, "implicit_age")
  }

  k <- k + rowSums(sapply(non_clim_vars, function(var) pars[[var]] * data[[var]] * implicit_age))

  for (clim_var in intersect(climatic_vars, pars_chosen)) {
    years <- seq(2019, 1985, by = -1)
    clim_columns <- paste0(clim_var, "_", years)
    k <- k + rowSums(sapply(clim_columns, function(col) pars[[clim_var]] * data[[col]]))
  }

  if ("B0" %in% pars) {
    return(pars[["B0"]] * (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
  } else {
    k[which(k > 7)] <- 7
    return(data[["mature_biomass"]] * (1 - exp(-k))^pars[["theta"]])
  }
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
      x = data$agbd - growth_curve(pars, data, pars_chosen),
      mean = 0,
      sd = pars["sd"],
      log = TRUE
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

  current_pars <- c(pars_basic, setNames(0.1, pars_chosen[1]))

  # Function to perform optimization
  optimize <- function(pars, pars_chosen_subset) {
    optim(pars,
      likelihood,
      fun = fun,
      data = data,
      pars_chosen = pars_chosen_subset,
      conditions = conditions
    )
  }

  o <- optimize(current_pars, pars_chosen[1])

  if (length(pars_chosen) > 1) {
    for (i in 2:length(pars_chosen)) {
      current_pars <- c(o$par, setNames(0.1, pars_chosen[i]))
      o <- optimize(current_pars, pars_chosen[1:i])
    }
  }

  pred <- growth_curve(o$par, data, pars_chosen)
  rsq <- cor(data$agbd, pred)^2
  print(paste("R-squared:", rsq))
  return(o)
}

################### Running model ###################

datafiles <- list(
  "data/5y_LULC.csv",
  "data/10y_LULC.csv",
  "data/15y_LULC.csv"
  # "data/all_LULC.csv"
)

dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

prop_amaz <- lapply(dataframes, function(df) {
  nrow(df %>% filter(biome == 1)) / nrow(df)
})
# 84% of them are in the Amazon anyways.

# adding names can make indexing complicated, so they are being referenced as a vector
names_dataframes <- c("data_5", "data_10", "data_15", "data_all")

if (run_all || run_one) {
  # Define conditions
  conditions <- list(
    'pars["theta"] > 10',
    'pars["theta"] < 0'
  )

  # in order to make them comparable, we only fit the columns that are present in all dataframes
  colnames_intersect <- Reduce(intersect, map(dataframes, colnames))
  # also exclude columns that are not predictors
  colnames_filtered <- colnames_intersect[!grepl(
    "agbd|latitude|longitude|prec|si|mature_biomass",
    colnames_intersect
  )]

  configurations <- list(
    c("age"),
    c("num_fires_before_regrowth"),
    c("age", "num_fires_before_regrowth"),
    c("num_fires_before_regrowth", "implicit_age", "k0"),
    c("indig", "protec", "cwd", "mean_prec", "mean_si", names(data)[str_detect(names(data), "LU")]),
    colnames_filtered,
    c(colnames_filtered, climatic_vars)
  )

  # adding names can make indexing complicated, so they are being referenced as a vector
  names_configurations <- c("age", "fires", "age_fires", "implicit_fire", "all_cat", "all_non_hist", "all_hist")
}

if (run_one) {
  data <- dataframes[[3]]
  pars_chosen <- configurations[[4]]

  # intercept, shape term, growth rate
  pars_basic <- c(theta = 5)

  if ("age" %in% pars_chosen) {
    conditions <- c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
  } else if (!("k0" %in% pars_chosen)) {
    pars_basic <- c(pars_basic, B0 = 40)
    conditions <- c(conditions, list('pars["B0"] < 0'))
}

  val <- run_optimization("nls", pars_basic, data, pars_chosen, conditions)
}


if (run_all) {
  # sum_squares_fit <- list()
  # pars_fit <- list()
  # Run optimization
  configurations <- configurations[1:4]
  dataframes <- dataframes[1:4]
  iterations <- expand.grid(seq_along(configurations), seq_along(dataframes))

  results <- foreach(i = iterations$Var1, j = iterations$Var2, .combine = "c") %dopar% {
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
    # pred <- growth_curve(o_iter$par, dataframes[[j]], configurations[[i]])
    # rsq <- rsq_fun(dataframes[[j]]$agbd, pred)
    # print(rsq)
    o_iter
  }

  # sum_squares_fit[[paste(names_dataframes[j], names_configurations[i])]] <- rsq
  # pars_fit[[paste(names_dataframes[j], names_configurations[i])]] <- o_iter$par
  # print(min(unlist(sum_squares_fit)))
}


run_growth_model <- function(data, initial_pars) {
  conditions <- list(
    'pars["theta"] > 10',
    'pars["theta"] < 0',
    'pars["B0"] < 0'
  )

  use_asymptote <- length(initial_pars) > 3

  if (use_asymptote) {
    conditions <- c(conditions, 'pars["B0"] > pars["A"]')
  }

  growth_curve <- if (use_asymptote) {
    function(pars, data) {
      pars[["B0"]] + (pars[["A"]] - pars[["B0"]]) * (1 - exp(-pars[["age"]] * data$age))^pars[["theta"]]
    }
  } else {
    function(pars, data) {
      pars[["B0"]] + (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-pars[["age"]]))^pars[["theta"]]
    }
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

# Usage
for (i in seq_along(dataframes)) {
  print("----------------------------------------------------")
  print(names_dataframes[i])

  # Without asymptote (using mature_biomass)
  result1 <- run_growth_model(dataframes[[i]], c(B0 = 40, theta = 5, age = 0.1))
  print(paste("R-squared (fixed asymptote, fit growth rate):", result1$r_squared))

  # # With asymptote
  # result2 <- run_growth_model(dataframes[[i]], c(B0 = 40, theta = 5, k = 0.1, A = 100))
  # print(paste("R-squared (fit asymptote, rate fit from age):", result2$r_squared))
}
