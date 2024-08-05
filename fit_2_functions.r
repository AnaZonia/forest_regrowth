################################################################################
#                                                                              #
#                 Forest Regrowth Model Functions and Utilities                #
#                                                                              #
#                            Ana Avila - August 2024                           #
#                                                                              #
#     This script defines the core functions used in the forest regrowth       #
#     modeling process (fit_3_run_model.r)                                     #
#                                                                              #
#     Functions included:                                                      #
#     - growth_curve                                                           #
#     - likelihood                                                             #
#     - filter_test_data                                                       #
#     - run_optim                                                              #
#     - run_gam                                                                #
#     - run_lm                                                                 #
#     - run_rf                                                                 #
#     - process_row                                                            #
#     - run_foreach                                                            #
#                                                                              #
################################################################################

library(mgcv)
library(randomForest)
set.seed(1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------ Growth Curve -------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------------------------- Likelihood -------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----------------- Ensuring testing data is within testing data range -------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

filter_test_data <- function(train_data, test_data) {
  # Calculate min and max for each column in train_data
  train_min <- sapply(train_data, min)
  train_max <- sapply(train_data, max)

  # Function to check if a row is within the range
  is_within_range <- function(row) {
    all(row >= train_min & row <= train_max)
  }

  # Apply the function to each row of test_data
  filtered_test_data <- test_data[apply(test_data, 1, is_within_range), ]

  return(filtered_test_data)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Optimization ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

  iter_pars <- c(
    theta = 5,
    setNames(
      rep(0.1, length(pars_basic)),
      c(pars_basic)
    ),
    setNames(
      rep(0, length(pars_chosen)),
      c(pars_chosen)
    )
  )

  conditions_iter <- conditions
  if ("age" %in% pars_basic) {
    conditions_iter <- c(conditions_iter, list('pars["age"] < 0', 'pars["age"] > 5'))
  }
  if ("B0" %in% pars_basic) {
    iter_pars["B0"] <- 40
    conditions_iter <- c(conditions_iter, list('pars["B0"] < 0'))
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FITTING
  # theta + other essential parameters
  remaining <- 1:length(iter_pars) # after find each one, remove from remaining list
  best <- list()
  best$AIC <- 0
  best$par <- c(iter_pars[c(1:(length(pars_basic) + 1))]) # identify which parameter was optimal.

  tmpinfo <- iter_pars
  tmpinfo[names(iter_pars)] <- NA
  tmpinfo <- c(likelihood = 0, tmpinfo)

  tmp <- as.data.frame(matrix(tmpinfo, nrow = 1))
  names(tmp) <- names(tmpinfo)
  taken <- length(remaining) + 1 # out of the range of values such that remaining[-taken] = remaining for the first iteration

  # iterate through # variables - keep best one, and go to next
  for (i in (length(pars_basic) + 2):length(iter_pars)) {

    tmp1 <- foreach(j = remaining[-taken]) %dopar% {
      inipar <- c(iter_pars[j], best$par) # as starting point, taking the best values from last time
      model <- optim(inipar,
        likelihood,
        fun = fun,
        data = train_data,
        conditions = conditions_iter
      )

      tmpinfo[names(inipar)] <- model$par
      tmpinfo["likelihood"] <- model$value
      return(tmpinfo) # standardize output to make easier to combine
    }

    tmp <- as.data.frame(do.call(rbind, tmp1))
    pos <- which.min(tmp$likelihood)
    tmp_AIC <- 2 * tmp$likelihood[pos] + 2 * i
    if (best$AIC == 0 | tmp_AIC < best$AIC) # keep the parameter values
      {
        best$AIC <- tmp_AIC
        # rewrite
        best$par <- tmp[pos, names(iter_pars)]
        taken <- which(!is.na(best$par))
        best$par <- best$par[taken]
        print(best$par)

      } else {
        print("No improvement. Exiting loop.")
        print(best$par)
        break
    }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  filtered_test_data <- filter_test_data(train_data, test_data)
  pred <- growth_curve(best$par, filtered_test_data)
  rsq <- cor(filtered_test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))
  
  return(list(
    pars = best$par,
    rsq = rsq
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------ Linear Model, Random Forest, and GAM -------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Output R-squared value and model results
  filtered_test_data <- filter_test_data(train_data, test_data)
  pred <- predict(model, newdata = filtered_test_data)
  rsq <- cor(filtered_test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))

  return(list(
    rsq = rsq
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Linear Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_lm <- function(train_data, test_data, pars_chosen) {
  lm_formula <- as.formula(paste("agbd ~", paste(pars_chosen, collapse = " + ")))

  model <- lm(lm_formula, data = train_data)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check for rank deficiency, and remove variables if necessary
  # (usually ones with too few unique values)
  aliased_vars <- summary(model)$aliased

  if (any(aliased_vars)) {
    problematic_vars <- names(aliased_vars)[aliased_vars]
    print(paste("Removing rank-deficient variables:", paste(problematic_vars, collapse = ", ")))

    # Remove problematic variables from the formula
    pars_chosen <- pars_chosen[!pars_chosen %in% problematic_vars]
    lm_formula <- as.formula(paste("agbd ~", paste(pars_chosen, collapse = " + ")))
    model <- lm(lm_formula, data = train_data)
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Output R-squared value and model results
  filtered_test_data <- filter_test_data(train_data, test_data)
  pred <- predict(model, newdata = filtered_test_data)
  rsq <- cor(filtered_test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))

  return(list(
    pars = t(summary(model)$coefficients[-1, 1, drop = FALSE]), # -1 to remove (Intercept),
    rsq = rsq
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random Forest
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_rf <- function(train_data, test_data, pars_chosen) {
  rf_formula <- as.formula(paste("agbd ~", paste(pars_chosen, collapse = " + ")))

  # reduce number of rows to 20,000 for faster computation
  train_data <- train_data[sample(nrow(train_data), 20000), ]

  model <- randomForest(rf_formula,
    data = train_data,
    ntree = 100, mtry = 2, importance = TRUE,
    keep.forest = TRUE, oob.score = TRUE, do.trace = 10, parallel = TRUE
  )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Output R-squared value and model results
  filtered_test_data <- filter_test_data(train_data, test_data)
  pred <- predict(model, newdata = filtered_test_data)
  rsq <- cor(filtered_test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))

  return(list(
    pars = t(importance(model)[, 1]),
    rsq = rsq
  ))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------------------- K-fold cross-validation ---------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cross_valid <- function(data, run_function, pars_chosen, pars_basic = NULL, conditions = NULL) {
  
  r_squared_values <- c()

  indices <- sample(c(1:5), nrow(data), replace = TRUE)

  for (i in 1:5) {
    # Define the test and train sets
    test_data <- data[indices == i, ]
    train_data <- data[!indices == i, ]

    # Run the function on the training data and evaluate on the test data
    if (identical(run_function, run_optim)) {
      model_output <- run_function("nls", pars_basic, pars_chosen, train_data, test_data, conditions)
    } else {
      model_output <- run_function(train_data, test_data, pars_chosen)
    }

    if (is.na(model_output$rsq)) {
      print(paste0("RSQ NA!", model_output$pars))
    } else {
      # Store the R squared value
      r_squared_values[i] <- model_output$rsq
    }

  }

  mean_r_squared <- mean(r_squared_values)
  sd_r_squared <- sd(r_squared_values)
  # return parameters of last fit
  return(list(rsq = mean_r_squared, rsq_sd = sd_r_squared, pars = model_output$pars))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------- Apply all models to the data ------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define a function to run foreach for different models
run_foreach <- function(iterations, model_type, run_function, conditions = NULL) {

  results <- foreach(iter = if (test_switch) test_rows else 1:nrow(iterations),
  .combine = "bind_rows", .packages = c("dplyr", "randomForest")) %dopar% {
    
    print(iter)
    i <- iterations$interval[iter]
    j <- iterations$data_par[iter]


    if (split_biome) {
      h <- iterations$biome[iter]
      biome_name <- biomes[[h]]
      data <- dataframes[[i]][[h]]
    } else {
      biome_name <- NULL
      data <- dataframes[[i]]
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

    if (model_type == "optim") {
      k <- iterations$basic_par[iter]
      pars_chosen <- data_pars[[j]]
      pars_basic <- basic_pars[[k]]
      pars_names <- data_pars_names[[j]]
      basic_pars_names <- basic_pars_names[[k]]
    } else {
      pars_chosen <- data_pars_lm[[j]]
      pars_names <- data_pars_lm_names[[j]]
      pars_basic <- NULL
      basic_pars_names <- NULL
    }
    
    cross_valid_output <- cross_valid(data, run_function, pars_chosen, pars_basic, conditions)

    # Organize result
    row <- process_row(cross_valid_output, model_type, intervals[[i]], pars_names,
      basic_pars_names = basic_pars_names, biome_name = biome_name
    )

    print(paste("Time so far: ", as.numeric(difftime(
      Sys.time(), start_time, units = "mins")), " minutes"))
    print(row)

    row
  
  }

  # Calculate and print the total time taken
  print(paste(
    model_type, "finished! Time for the whole operation: ",
    as.numeric(difftime(Sys.time(), start_time, units = "hours")), " hours"
  ))

  # Combine all results into a single dataframe
  df <- as.data.frame(results)

  # Write the dataframe to a CSV file
  # if (!all_switch) {
  if (split_biome) {
    write.csv(df, paste0("./data/", region, "_results_", model_type, "_split_biome.csv"), row.names = FALSE)
  } else {
    write.csv(df, paste0("./data/", region, "_results_", model_type, ".csv"), row.names = FALSE)
  }
  # }
  
  return(df)

}    


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#--------------------- Prepare iteration results as dataframe row -----------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Intakes:
#   output <- the output of the model as a matrix (coefficients or rf importance)
#   data_name <- "5y", "10y", "15y", or "all"
#   model_type <- "optim", "lm", "rf", or "gam"
#   rsq <- the calculated r squared for each iteration of the model

# Define helper functions
process_row <- function(output, model_type, data_name, data_pars_names, basic_pars_names = NULL,
biome_name = NULL) {

  row <- as.data.frame(output$pars)
  # Set up columns of parameters that are not included in this iteration as NA
  missing_cols <- setdiff(unique(unlist(c(data_pars_lm, non_data_pars, data_pars))), names(row))
  row[missing_cols] <- NA
  row <- row[, unique(unlist(c(data_pars_lm, non_data_pars, data_pars)))]
  
  row$data_name <- data_name
  row$model_type <- model_type
  row$data_pars <- data_pars_names
  row$rsq <- output$rsq
  row$rsq_sd <- output$rsq_sd

  if (model_type == "optim"){
    row$basic_pars <- basic_pars_names
  } else {
    row$basic_pars <- NA
  }

  # Define the desired order of columns
  desired_column_order <- c(
    "data_name", "data_pars", "model_type", "rsq", "rsq_sd", "basic_pars",
    "mature_biomass", "age"
  )

  if (!is.null(biome_name)) {
    row <- row %>%
      mutate(biome = biome_name) %>%
      select(biome, all_of(desired_column_order), all_of(non_data_pars), everything())
  } else {
    row <- row %>%
      select(all_of(desired_column_order), all_of(non_data_pars), everything())
  }

  return(row)
}
