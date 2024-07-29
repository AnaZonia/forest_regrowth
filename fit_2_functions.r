################################################################################
#                                                                              #
#                 Forest Regrowth Model Functions and Utilities                #
#                                                                              #
#                              Ana Avila - July 2024                           #
#                                                                              #
#     This script defines the core functions used in the forest regrowth       #
#     modeling process (fit_3_run_model.r)                                     #
#                                                                              #
#     Functions included:                                                      #
#     - growth_curve                                                           #
#     - likelihood                                                             #
#     - run_optimization                                                       #
#     - run_gam_lm                                                             #
#     - run_rf                                                                 #
#     - process_row                                                            #
#     - run_foreach                                                            #
#                                                                              #
################################################################################

library(mgcv)
library(randomForest)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------------- Growth Curve -----------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#     This function calculates the Chapman-Richards growth curve based on the provided
#     parameters and data. It can incorporate:
#     - Yearly-changing climatic parameters (prec and si)
#     - The intercept term defined either as B0 or B0_exp (out or in of the exponential)
#     - A growth rate intercept term k0
#     - Ages of secondary forests as an explicit parameter (part of "pars")
#       or as an implicit parameter (multiplying all non-yearly predictors by age)
#
# Intakes:
#   pars <- a vector with the named parameters to be included
#   data <- the chosen training dataframe

growth_curve <- function(pars, data) {

  # Define parameters that are not expected to change yearly (not prec or si)
  non_clim_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars))
  # Define whether age is an explicit or implicit parameter (to multiply the other parameters by)
  implicit_age <- if (!"age" %in% names(pars)) data[["age"]] else rep(1, nrow(data))
  # Define whether the intercept k0 is to be included in the growth rate k
  k <- if ("k0" %in% names(pars)) pars[["k0"]] * implicit_age else rep(0, nrow(data))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate the growth rate k
  k <- k + rowSums(sapply(non_clim_pars, function(par) pars[[par]] * data[[par]] * implicit_age))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Add yearly-changing climatic parameters to the growth rate k (if included in the parameter set)
  for (clim_par in intersect(climatic_pars, names(pars))) {
    years <- seq(2019, 1985, by = -1)
    clim_columns <- paste0(clim_par, "_", years)
    k <- k + rowSums(sapply(clim_columns, function(col) pars[[clim_par]] * data[[col]]))
  }

  k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)

  if ("B0" %in% names(pars)) {
    return(pars[["B0"]] * (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
  } else {
    return(data[["mature_biomass"]] * (1 - exp(-(pars[["B0_exp"]] + k)))^pars[["theta"]])
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------------- Likelihood -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Intakes:
#   fun <- the function to be used, either "nls" (Nonlinear Least Squares)
#          or "nll" (Negative Log Likelihood)
#   pars <- a vector with the named parameters to be included
#   data <- the chosen training dataframe
#   conditions <- ranges of parameters to be restricted to

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

  # Check whether any of the parameters is breaking the conditions (e.g. negative values)
  if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
    return(-Inf)
  } else if (is.na(result) || result == 0) {
    return(-Inf)
  } else {
    return(result)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----------------- Optimization -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#     This function prepares the optimization process for the forest regrowth model.
#     It runs the optimization process iteratively, incorporating one parameter of the
#     chosen set at a time.
#
#     It also determines:
#     - A new "sd" parameter if the function is set to "nll" (Negative Log Likelihood)
#     - The theta ("shape") parameter of the Chapman-Richards growth curve
#     - The values of the initial parameters (0.1 for all except for theta and B0)
#
# Intakes:
#   fun <- the function to be used, either "nls" (Nonlinear Least Squares)
#          or "nll" (Negative Log Likelihood)
#   pars_basic <- a vector with the named parameters that do not correspond to data,
#                 but are to be included as part of the functional form
#   pars_chosen <- a vector with the chosen parameter set to be included
#                  (corresponding to columns in dataframe)
#   train_data <- the chosen training dataframe
#   test_data <- the chosen testing dataframe (for rsq calculation)
#   conditions <- ranges of parameters to be restricted to

run_optim <- function(fun, pars_basic, pars_chosen, train_data, test_data, conditions) {
  # Define function to perform optimization with input parameters
  optimize <- function(pars, data) {
    optim(pars,
      likelihood,
      fun = fun,
      data = data,
      conditions = conditions
    )
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define initial parameters for the first iteration
  first_iter_pars <- c(
    theta = 5,
    setNames(
      rep(0.1, length(c(pars_basic, pars_chosen[1]))),
      c(pars_basic, pars_chosen[1])
    )
  )

  if (fun == "nll") {
    first_iter_pars <- c(first_iter_pars, setNames(0.1, "sd"))
  }

  if ("B0" %in% pars_basic) {
    first_iter_pars["B0"] <- 40
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  optimize_iteratively <- function(data) {
    # Run first iteration of optimization
    model <- optimize(first_iter_pars, data)

    # Intake parameters of previous iteration and run optimization for the next parameter
    if (length(pars_chosen) > 1) {
      for (i in 2:length(pars_chosen)) {
        current_pars <- c(model$par, setNames(0.1, pars_chosen[i]))
        model <- optimize(current_pars, data)
      }
    }

    return(model)
  }

  model <- optimize_iteratively(train_data)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Output R-squared value and model results
  pred <- growth_curve(model$par, test_data)
  rsq <- cor(test_data$agbd, pred)^2

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # fitting optim with testing and training data raised issues.
  # A few times, the predictions for the training data were fine, but for the testing data there were NAs.
  # In those cases, I told optim to run again with the entire dataframe and let me know if there was a discrepancy.
  discrepancy = FALSE
  if (is.na(rsq)){
    discrepancy = TRUE
    all_data <- bind_rows(train_data, test_data)
    model <- optimize_iteratively(all_data)
    pred <- growth_curve(model$par, all_data)
    rsq <- cor(all_data$agbd, pred)^2
  }

  print(paste("R-squared:", rsq))
  return(list(
    model = model,
    rsq = rsq,
    discrepancy = discrepancy
  ))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------- Linear Model, Random Forest, and GAM -----------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Intake:
#   train_data <- the chosen training dataframe
#   test_data <- the chosen testing dataframe (for rsq calculation)
#   pars_chosen <- a vector with the chosen parameter set to be included
#                  (corresponding to columns in dataframe)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GAM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_gam <- function(train_data, test_data, pars_chosen) {
  # Separate pars_chosen into continuous and categorical variables
  continuous <- c(pars_chosen[!pars_chosen %in% pars_categ])
  categorical <- pars_chosen[pars_chosen %in% pars_categ]

  formula <- as.formula(paste(
    "agbd ~",
    paste(sapply(continuous, function(x) paste0("s(", x, ")")), collapse = " + "),
    "+",
    paste(categorical, collapse = " + ")
  ))

  model <- gam(formula, data = train_data)
  pred <- predict(model, newdata = test_data)
  rsq <- cor(test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))
  return(list(
    model = model,
    rsq = rsq
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Linear Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_lm <- function(train_data, test_data, pars_chosen) {
  rf_lm_formula <- as.formula(paste("agbd ~", paste(pars_chosen, collapse = " + ")))

  model <- lm(rf_lm_formula, data = train_data)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check for rank deficiency, and remove variables if necessary
  # (usually ones with too few unique values)
  aliased_vars <- summary(model)$aliased

  if (any(aliased_vars)) {
    problematic_vars <- names(aliased_vars)[aliased_vars]
    print(paste("Removing rank-deficient variables:", paste(problematic_vars, collapse = ", ")))

    # Remove problematic variables from the formula
    pars_chosen <- pars_chosen[!pars_chosen %in% problematic_vars]
    rf_lm_formula <- as.formula(paste("agbd ~", paste(pars_chosen, collapse = " + ")))
    model <- lm(rf_lm_formula, data = train_data)
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Output R-squared value and model results
  pred <- predict(model, newdata = test_data)
  rsq <- cor(test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))
  return(list(
    model = model,
    rsq = rsq
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random Forest
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_rf <- function(train_data, test_data, pars_chosen) {
  rf_lm_formula <- as.formula(paste("agbd ~", paste(pars_chosen, collapse = " + ")))

  # reduce number of rows to 20,000 for faster computation
  train_data <- train_data[sample(nrow(train_data), 20000), ]

  model <- randomForest(rf_lm_formula,
    data = train_data,
    ntree = 100, mtry = 2, importance = TRUE,
    keep.forest = TRUE, oob.score = TRUE, do.trace = 10, parallel = TRUE
  )

  pred <- predict(model, newdata = test_data)
  rsq <- cor(test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))
  return(list(
    model = importance(model),
    rsq = rsq
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------- Prepare iteration results as dataframe row -----------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define a function to run foreach for different models
run_foreach <- function(iterations, model_type, run_function, conditions = NULL) {

  results <- foreach(iter = if (test_switch) 1 else 1:nrow(iterations),
  .combine = "bind_rows", .packages = c("dplyr", "randomForest")) %dopar% {
    
    i <- iterations$dataframe[iter]
    j <- iterations$data_par[iter]

    train_data <- train_dataframes[[i]]
    test_data <- test_dataframes[[i]]

    if (model_type == "optim") {
      k <- iterations$basic_par[iter]
      pars_chosen <- data_pars[[j]]
      pars_basic <- basic_pars[[k]]
      pars_names <- data_pars_names[[j]]
      basic_pars_names <- basic_pars_names[[k]]

      conditions_iter <- conditions
      if ("age" %in% pars_chosen) {
        conditions_iter <- c(conditions_iter, list('pars["age"] < 0', 'pars["age"] > 5'))
      } else if ("B0" %in% pars_chosen) {
        conditions_iter <- c(conditions_iter, list('pars["B0"] < 0'))
      }


      # Call run_optim with the correct parameters
      model_output <- run_function("nls", pars_basic, pars_chosen, train_data, test_data, conditions)
    } else {
      pars_chosen <- data_pars_lm[[j]]
      pars_names <- data_pars_lm_names[[j]]
      basic_pars_names <- NULL

      # Call run_function for non-optim models
      model_output <- run_function(train_data, test_data, pars_chosen)
    }

    # Extract fitting parameters based on model type
    fit_pars <- switch(model_type,
      "optim" = model_output$model$par,
      "lm" = summary(model_output$model)$coefficients[-1, 1, drop = FALSE],
      "rf" = model_output$model[, 1]
    )

    # Organize result
    row <- process_row(fit_pars, intervals[[i]], model_type, pars_names, model_output,
      basic_pars_names = basic_pars_names
    )

    print(paste("Time so far: ", as.numeric(difftime(Sys.time(), start_time, units = "mins")), " minutes"))
    
    row
  
  }

  # Calculate and print the total time taken
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
  print(paste(
    model_type, "finished! Time for the whole operation: ",
    total_time, " hours"
  ))

  # Combine all results into a single dataframe
  df <- as.data.frame(results)

  # Write the dataframe to a CSV file
  write.csv(df, paste0("./data/results_", model_type, ".csv"), row.names = FALSE)

  return(df)

}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------- Prepare iteration results as dataframe row -----------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Intakes:
#   output <- the output of the model as a matrix (coefficients or rf importance)
#   data_name <- "5y", "10y", "15y", or "all"
#   model_type <- "optim", "lm", "rf", or "gam"
#   rsq <- the calculated r squared for each iteration of the model

# Define helper functions
process_row <- function(fit_pars, data_name, model_type, data_pars_names, output, basic_pars_names = NULL) {
  row <- as.data.frame(t(fit_pars))
  # Set up columns of parameters that are not included in this iteration as NA
  missing_cols <- setdiff(unique(unlist(c(data_pars_lm, non_data_pars, data_pars))), names(row))
  row[missing_cols] <- NA
  row <- row[, unique(unlist(c(data_pars_lm, non_data_pars, data_pars)))]
  
  row$data_name <- data_name
  row$model_type <- model_type
  row$data_pars <- data_pars_names
  row$rsq <- output$rsq

  if (model_type == "optim"){
    row$discrepancy <- output$discrepancy
    row$basic_pars <- basic_pars_names
  } else {
    row$discrepancy <- NA
    row$basic_pars <- NA
  }
  
  row %>% select(data_name, "data_pars", "basic_pars", model_type, rsq, "discrepancy", "mature_biomass", "age", all_of(non_data_pars), everything())
}

