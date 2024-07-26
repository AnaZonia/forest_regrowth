####################################################################
############### Process data and define functions ##################
####################################################################
# Ana Avila - July 2024
# ~~~~~~~~~~~~~~~~~~~~
# Intakes:
# Outputs:
####################################################################
library(ggplot2)
library(foreach)
library(doParallel)

source("fit_1_import_data.r")

#------------------ SWITCHES ------------------#

run_all <- TRUE
run_one <- TRUE
fit_gaus_ker <- FALSE

#------------------ FUNCTIONS ------------------#


# Intakes pars <- a vector with the initial parameters to be included
# data <- the dataframe with the predictors
# pars_chosen <- the list of parameters to be added into the shape term


growth_curve <- function(pars, data, pars_chosen) {

  implicit_age <- 1
  if ("k0" %in% pars_chosen) {
    implicit_age <- data[["age"]]
    k <- pars[["k0"]] * implicit_age
    non_clim_pars <- setdiff(pars_chosen, c("k0", "implicit_age", climatic_pars))
  } else {
    k <- rep(0, nrow(data))
    non_clim_pars <- setdiff(pars_chosen, climatic_pars)
  }

  k <- k + rowSums(sapply(non_clim_pars, function(par) pars[[par]] * data[[par]] * implicit_age))
  
  for (clim_par in intersect(climatic_pars, pars_chosen)) {
    years <- seq(2019, 1985, by = -1)
    clim_columns <- paste0(clim_par, "_", years)
    k <- k + rowSums(sapply(clim_columns, function(col) pars[[clim_par]] * data[[col]]))
  }

  if ("B0" %in% pars) {
    return(pars[["B0"]] * (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
  } else {
    k[which(k > 7)] <- 7
    return(data[["mature_biomass"]] * (1 - exp(-(pars[["B0_exp"]] + k)))^pars[["theta"]])
  }
}

growth_curve(pars, data, pars_chosen)
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

#------------------ Global Variables ------------------#

climatic_pars <- c("prec", "si")

intervals <- list(
  "5y",
  "10y",
  "15y",
  "all"
)

datafiles <- lapply(intervals, function(file) {
  paste0("./data/", file, "_LULC_dist_amaz_500.csv")
})
dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

# adding names can make indexing complicated, so they are being referenced as a vector
names_dataframes <- c("data_5", "data_10", "data_15", "data_all")

################### Running model ###################


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
    "age|agbd|latitude|longitude|prec|si|mature_biomass|distance|biome|cwd",
    colnames_intersect
  )]

  configurations <- list(
    c(),
    c("cwd", "mean_prec", "mean_si"),
    c("cwd", climatic_pars),
    colnames_filtered[!grepl("ecoreg|soil", colnames_filtered)],
    colnames_filtered,
    c(colnames_filtered, "cwd", "mean_prec", "mean_si"),
    c(colnames_filtered, climatic_pars),
  )

  basic_pars <- list(
    c(),
    c("implicit_age", "k0", "B0_exp"),
    c("implicit_age", "k0", "B0"),
    c("age", "B0")
  )

  combined_list <- list()
  for (config in configurations) {
    for (basic in basic_pars) {
      combined_list <- c(combined_list, list(c(config, basic)))
    }
  }

  # adding names can make indexing complicated, so they are being referenced as a vector
  names_configurations <- c(
    "none", "clim_mean", "clim_hist", "lulc",
    "lulc_ecoreg_soil", "lulc_ecoreg_soil_clim_mean", "lulc_ecoreg_soil_clim_hist"
  )
}

run_one = TRUE
if (run_one) {

  pars_basic <- c(theta = 5)

  # intercept, shape term, growth rate
  data <- dataframes[[3]]

  if ("age" %in% pars_chosen) {
    conditions_iter <- c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
  } else if (!("k0" %in% pars_chosen)) {
    pars_basic_iter <- c(pars_basic, B0 = 40)
    conditions_iter <- c(conditions, list('pars["B0"] < 0'))
  }
  
  val <- run_optimization("nls", pars_basic, data, pars_chosen, conditions)
}
  pars_basic <- c(theta = 5, k0 = 0.1, B0_exp = 0.01)


if (run_all) {
  # sum_squares_fit <- list()
  # pars_fit <- list()
  # Run optimization
  configurations <- configurations[1:length(configurations)]
  dataframes <- dataframes[1:length(dataframes)]
  iterations <- expand.grid(seq_along(configurations), seq_along(dataframes))
  registerDoParallel(cores = 15)

  pars_basic <- c(theta = 5, k0 = 0.1, B0_exp = 0.01)

  results <- foreach(i = iterations$Var1, j = iterations$Var2, .combine = "c") %dopar% {
    pars_basic_iter <- pars_basic
    conditions_iter <- conditions
    if ("age" %in% configurations[[i]]) {
      conditions_iter <- c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
    } else if ("B0" %in% pars_chosen) {
      # pars_basic_iter <- c(pars_basic, B0 = 40)
      conditions_iter <- c(conditions, list('pars["B0"] < 0'))
    }

    o_iter <- run_optimization(
      "nls", pars_basic_iter, dataframes[[j]], configurations[[i]], conditions_iter
    )

    print(names_dataframes[j])
    print(names_configurations[i])
    o_iter
  }

  saveRDS(results, file = "results.rds")
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

tst <- readRDS("results.rds")
tst
# # Usage
# for (i in seq_along(dataframes)) {
#   print("----------------------------------------------------")
#   print(names_dataframes[i])

#   # Without asymptote (using mature_biomass)
#   result1 <- run_growth_model(dataframes[[i]], c(B0 = 40, theta = 5, age = 0.1))
#   print(paste("R-squared (fixed asymptote, fit growth rate):", result1$r_squared))

#   # # With asymptote
#   # result2 <- run_growth_model(dataframes[[i]], c(B0 = 40, theta = 5, k = 0.1, A = 100))
#   # print(paste("R-squared (fit asymptote, rate fit from age):", result2$r_squared))
# }
