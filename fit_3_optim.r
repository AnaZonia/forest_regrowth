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
# pars <- the list of parameters to be added into the shape term

growth_curve <- function(pars, data) {

  non_clim_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars))

  implicit_age <- if ("implicit_age" %in% names(pars)) data[["age"]] else rep(1, nrow(data))
  
  k <- if ("k0" %in% names(pars)) pars[["k0"]] * implicit_age else rep(0, nrow(data))

  k <- k + rowSums(sapply(non_clim_pars, function(par) pars[[par]] * data[[par]] * implicit_age))

  for (clim_par in intersect(climatic_pars, names(pars))) {
    years <- seq(2019, 1985, by = -1)
    clim_columns <- paste0(clim_par, "_", years)
    k <- k + rowSums(sapply(clim_columns, function(col) pars[[clim_par]] * data[[col]]))
  }
  
  k[which(k > 7)] <- 7 # in order to avoid increasinly small values for exp(k)

  if ("B0" %in% names(pars)) {
    return(pars[["B0"]] * (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
  } else {
    return(data[["mature_biomass"]] * (1 - exp(-(pars[["B0_exp"]] + k)))^pars[["theta"]])
  }
}

# Calculates Nonlinear Least Squares
# Intakes:
# fun <- the function to be used, either "nls" or "nll"
# pars <- a vector with the initial parameters to be included
# data <- the dataframe with the predictors
# pars_chosen <- the list of parameters to be added into the shape term
# conditions <- ranges of parameters to be restricted to

likelihood <- function(fun, pars, data, conditions) {
  if (fun == "nll") {
    result <- -sum(dnorm(
      x = data$agbd - growth_curve(pars, data),
      mean = 0,
      sd = pars["sd"],
      log = TRUE
    ), na.rm = TRUE)
    conditions <- c(conditions, 'pars["sd"] < 0')
  } else {
    result <- sum((growth_curve(pars, data) - data$agbd)^2)
  }

  if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
    return(-Inf)
  } else if (is.na(result) || result == 0) {
    return(-Inf)
  } else {
    return(result)
  }
}

run_optimization <- function(fun, pars_basic, data_pars, data, conditions) {

  # Run optimization
  if (fun == "nll") {
    pars <- c(pars, setNames(0.1, "sd"))
  }

  current_pars <- c(
    theta = 5,
    setNames(
      rep(0.1, length(c(pars_basic, data_pars[1]))),
      c(pars_basic, data_pars[1])
    )
  )

  if ("B0" %in% pars_basic) {
    current_pars["B0"] <- 40
  }
  
  # Function to perform optimization
  optimize <- function(pars) {
    optim(pars,
      likelihood,
      fun = fun,
      data = data,
      conditions = conditions
    )
  }

  o <- optimize(current_pars)

  if (length(data_pars) > 1) {
    for (i in 2:length(data_pars)) {
      current_pars <- c(o$par, setNames(0.1, data_pars[i]))
      o <- optimize(current_pars)
    }
  }

  pred <- growth_curve(o$par, data)
  rsq <- cor(data$agbd, pred)^2
  print(paste("R-squared:", rsq))
  return(results <- list(
    model = o,
    rsq = rsq
  ))
}

#------------------ Global Variables ------------------#

climatic_pars <- c("prec", "si")
# parameters not fit with data - just for the functional form
non_data_pars <- c("implicit_age", "k0", "B0_exp", "B0", "theta")

#------------------ Import Data -------------------#

intervals <- list(
  "5y",
  "10y",
  "15y",
  "all"
)

datafiles <- lapply(intervals, function(file) {
  paste0("./data/", file, "_LULC_mat_dist.csv")
})

dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

dataframes <- lapply(dataframes, function(df) {
  sample_n(df, size = 60000)
})

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

  data_pars <- list(
    c("cwd", "mean_prec", "mean_si"),
    c("cwd", climatic_pars),
    colnames_filtered[!grepl("ecoreg|soil", colnames_filtered)],
    colnames_filtered,
    c(colnames_filtered, "cwd", "mean_prec", "mean_si"),
    c(colnames_filtered, climatic_pars)
  )

  basic_pars <- list(
    c("B0"),
    c("B0_exp"),
    c("implicit_age", "k0", "B0_exp"),
    c("implicit_age", "k0", "B0"),
    c("age", "k0", "B0"),
    c("age", "B0")
  )

  # adding names can make indexing complicated, so they are being referenced as a vector
  names_configurations <- c(
    "none", "clim_mean", "clim_hist", "lulc",
    "lulc_ecoreg_soil", "lulc_ecoreg_soil_clim_mean", "lulc_ecoreg_soil_clim_hist"
  )
}

if (run_one) {
  pars_chosen <- data_pars[[2]]
  pars_basic <- basic_pars[[1]]

  if ("age" %in% pars_chosen) {
    conditions_iter <- c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
  } else if ("B0" %in% pars_basic) {
    conditions_iter <- c(conditions, list('pars["B0"] < 0'))
  }
  
  o_iter <- run_optimization("nls", pars_basic, pars_chosen, data, conditions)

}


if (run_all) {
  # sum_squares_fit <- list()
  # pars_fit <- list()
  # Run optimization
  iterations <- expand.grid(
    dataframe = seq_along(dataframes),
    data_par = seq_along(data_pars),
    basic_par = seq_along(basic_pars)
  )

  registerDoParallel(cores = 15)
  
  results <- foreach(iter = 1:nrow(iterations),  .combine = 'rbind', .packages = c('dplyr')) %dopar% {
    
    i <- iterations$dataframe[iter]
    j <- iterations$data_par[iter]
    k <- iterations$basic_par[iter]

    data <- dataframes[[i]]
    pars_chosen <- data_pars[[j]]
    pars_basic <- basic_pars[[k]]
    
    if ("age" %in% pars_chosen) {
      conditions_iter <- c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
    } else if ("B0" %in% pars_chosen) {
      conditions_iter <- c(conditions, list('pars["B0"] < 0'))
    }

    o_iter <- run_optimization("nls", pars_basic, pars_chosen, data, conditions_iter)

    new_row <- o_iter$model$par
    new_row <- as.data.frame(t(new_row))

    new_row$model_name <- intervals[[i]]
    new_row$model_type <- "optim"
    new_row$rsq <- o_iter$rsq

    # Reorder columns to have model_name first and rsq second
    new_row <- new_row %>% select(model_name, model_type, rsq, everything())
    print(new_row)
    new_row
  }

  # Combine all results into a single dataframe
  lm_df <- as.data.frame(results)

  # Write the dataframe to a CSV file
  write.csv(lm_df, "./data/fit_results.csv", row.names = FALSE)

}

