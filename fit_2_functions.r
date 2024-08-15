# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              
#                 Forest Regrowth Model Functions and Utilities                
#                                                                              
#                            Ana Avila - August 2024                           
#                                                                              
#     This script defines the core functions used in the forest regrowth       
#     modeling process (fit_3_run_model.r)                                     
#                                                                              
#     Functions included:                                                      
#     - growth_curve                                                           
#     - likelihood                                                             
#     - filter_test_data                                                       
#     - run_optim                                                              
#     - run_gam                                                                
#     - run_lm                                                                 
#     - run_rf                                                                 
#     - process_row                                                            
#     - run_foreach                                                            
#                                                                              
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
  if (length(non_clim_pars) > 0){
    k <- k + rowSums(sapply(non_clim_pars, function(par) pars[[par]] * data[[par]] * implicit_age))
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Add yearly-changing climatic parameters to the growth rate k (if included in the parameter set)
  for (clim_par in intersect(climatic_pars, names(pars))) {
    years <- seq(2019, 1985, by = -1)
    clim_columns <- paste0(clim_par, "_", years)
    k <- k + rowSums(sapply(clim_columns, function(col) pars[[clim_par]] * data[[col]]))
  }

  # na_indices <- which(k < 0)
  # if (length(na_indices) > 0){
  #   # cat("number of rows removed due to negative k:", length(na_indices), "out of", length(k))
  #   k[na_indices] <- 0
  # }

  k[which(k > 7)] <- 7 # Constrains k to avoid increasinly small values for exp(k) (local minima at high k)
  
  if (fit_logistic) {
    return(data[["mature_biomass"]] * (1 / (1 + exp(k))))
    # return(pars[["B0"]] * data[["mature_biomass"]] * exp(k)) / ((data[["mature_biomass"]] - B0) + B0 * exp(k))
  } else {
    if ("B0" %in% names(pars)) {
      return(pars[["B0"]] + (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
    } else {
      return(data[["mature_biomass"]] * (1 - exp(-(pars[["B0_exp"]] + k)))^pars[["theta"]])
    }
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

rsq_filtered_data <- function(train_data, test_data, model, model_type) {
  # Calculate min and max for each column in train_data
  train_min <- sapply(train_data, min)
  train_max <- sapply(train_data, max)

  # Function to check if a row is within the range
  is_within_range <- function(row) {
    all(row >= train_min & row <= train_max)
  }

  # Apply the function to each row of test_data
  filtered_test_data <- test_data[apply(test_data, 1, is_within_range), ]

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Output R-squared value and model results
  if (model_type == "optim"){
    pred <- growth_curve(model$par, filtered_test_data)
  } else {
    pred <- predict(model, newdata = filtered_test_data)
  }

  rsq <- cor(filtered_test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))

  return(rsq)
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
#   basic_pars_iter <- a vector with the named parameters that do not correspond to data,
#                 but are to be included as part of the functional form
#   data_pars_iter <- a vector with the chosen parameter set to be included
#                  (corresponding to columns in dataframe)
#   train_data <- the chosen training dataframe
#   test_data <- the chosen testing dataframe (for rsq calculation)
#   conditions <- ranges of parameters to be restricted to


run_optim <- function(fun, pars, train_data, conditions, test_data = NULL) {

  if ("age" %in% names(pars)) {
    conditions <- c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
  }
  if ("B0" %in% names(pars)) {
    conditions <- c(conditions, list('pars["B0"] < 0'))
  }

  
  # Use tryCatch to handle errors during optimization
  model <- tryCatch(
    {
      optim(pars,
        likelihood,
        fun = fun,
        data = train_data,
        conditions = conditions
      )
    },
    error = function(e) {
      # Print the parameters if an error occurs
      print(e)
      print(pars)
    }
  )


  if (is.null(test_data)){
    return(model)
  } else {
    return(list(
      pars = t(model$par),
      rsq = rsq_filtered_data(train_data, test_data, model, "optim")
    ))
  }

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------ Linear Model, Random Forest, and GAM -------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Intake:
#   train_data <- the chosen training dataframe
#   test_data <- the chosen testing dataframe (for rsq calculation)
#   pars <- a vector with the chosen parameter set to be included
#                  (corresponding to columns in dataframe)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GAM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_gam <- function(train_data, test_data, pars) {

  # Define categorical parameters
  pars_categ <- c("indig", "protec", names(train_data)[str_detect(names(train_data), "LU|ecoreg|soil")])

  # Separate pars into continuous and categorical variables
  continuous <- pars[!pars %in% pars_categ]
  categorical <- pars[pars %in% pars_categ]

  formula <- as.formula(paste(
    "agbd ~",
    paste(sapply(continuous, function(x) paste0("s(", x, ")")), collapse = " + "),
    "+",
    paste(categorical, collapse = " + ")
  ))

  model <- gam(formula, data = train_data)

  return(list(
    rsq = rsq_filtered_data(train_data, test_data, model, "gam")
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Linear Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_lm <- function(train_data, test_data, pars) {
  lm_formula <- as.formula(paste("agbd ~", paste(pars, collapse = " + ")))

  model <- lm(lm_formula, data = train_data)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check for rank deficiency, and remove variables if necessary
  # (usually ones with too few unique values)
  aliased_vars <- summary(model)$aliased

  if (any(aliased_vars)) {
    problematic_vars <- names(aliased_vars)[aliased_vars]
    print(paste("Removing rank-deficient variables:", paste(problematic_vars, collapse = ", ")))

    # Remove problematic variables from the formula
    pars <- pars[!pars %in% problematic_vars]
    lm_formula <- as.formula(paste("agbd ~", paste(pars, collapse = " + ")))
    model <- lm(lm_formula, data = train_data)
  }


  return(list(
    pars = t(summary(model)$coefficients[-1, 1, drop = FALSE]), # -1 to remove (Intercept),
    rsq = rsq_filtered_data(train_data, test_data, model, "lm")
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random Forest
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_rf <- function(train_data, test_data, pars) {

  rf_formula <- as.formula(paste("agbd ~", paste(pars, collapse = " + ")))

  model <- randomForest(rf_formula,
    data = train_data,
    ntree = 100, mtry = 2, importance = TRUE,
    keep.forest = TRUE, oob.score = TRUE, do.trace = 10, parallel = TRUE
  )

  return(list(
    pars = t(importance(model)[, 1]),
    rsq = rsq_filtered_data(train_data, test_data, model, "rf")
  ))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------- Apply all models to the data ------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Intake:
#   train_data <- the chosen training dataframe
#   test_data <- the chosen testing dataframe (for rsq calculation)
#   data_pars_iter <- a vector with the chosen parameter set to be included
#                  (corresponding to columns in dataframe)


find_combination_pars <- function(iterations, fun = "nls") {

  ideal_par_combination <- list()

  for (iter in 1:nrow(iterations)) {

    i <- iterations$interval[iter]
    j <- iterations$data_par[iter]
    k <- iterations$basic_par[iter]

    data_pars_iter <- data_pars[[j]]
    basic_pars_iter <- basic_pars[[k]]

    if (split_biome) {
      h <- iterations$biome[iter]
      data <- dataframes[[i]][[h]]
    } else {
      data <- dataframes[[i]]
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
    all_pars_iter <- c(setNames(
      rep(0, length(data_pars_iter)),
      c(data_pars_iter)
    ))
    
    if (!fit_logistic) {
      all_pars_iter[["theta"]] <- 1
    }

    if ("B0" %in% basic_pars_iter) {
      all_pars_iter["B0"] <- mean(data[["agbd"]])
    }
    
    if ("age" %in% basic_pars_iter) {
      all_pars_iter["age"] <- 0
    }

    if ("k0" %in% basic_pars_iter) {
      if ("B0" %in% basic_pars_iter) {
        all_pars_iter["k0"] <- -log(1 - mean(data[["agbd"]]) / mean(data[["mature_biomass"]]))
      }
      all_pars_iter["B0_exp"] <- -log(1 - mean(data[["agbd"]]) / mean(data[["mature_biomass"]]))
      all_pars_iter["k0"] <- 0
    }


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Account for dummy variables
    # after find each one, remove from remaining list

    categorical <- c("ecoreg", "soil", "last_LU")
    
    for (cat_var in categorical) {
      dummy_indices <- grep(cat_var, data_pars_iter)
      if (length(dummy_indices) > 0) {
        data_pars_iter <- c(data_pars_iter[-dummy_indices], cat_var)
      }
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    remaining <- 1:length(data_pars_iter)
    best <- list()
    best$AIC <- 0
    # initial vector is theta + all basic parameters
    best$par <- all_pars_iter[names(all_pars_iter) %in% c(basic_pars_iter, "theta")]
    taken <- length(remaining) + 1 # out of the range of values such that remaining[-taken] = remaining for the first iteration
    
    base_row <- all_pars_iter
    base_row[names(all_pars_iter)] <- NA
    base_row <- c(likelihood = 0, base_row)

    # iterate through # variables - keep best one, and go to next
    for (i in 1:length(data_pars_iter)) {
      
      optim_remaining_pars <- foreach(j = remaining[-taken]) %dopar% {
        # check for categorical variables (to be included as a group)
        if (data_pars_iter[j] %in% c("last_LU", "ecoreg", "soil")) {
          all_pars_iter_var <- all_pars_iter[grep(data_pars_iter[j], names(all_pars_iter))]
          inipar <- c(best$par, all_pars_iter_var)
        } else {
          inipar <- c(best$par, all_pars_iter[data_pars_iter[j]]) # as starting point, taking the best values from last time
        }

        model <- run_optim("nls", inipar, data, conditions)

        iter_row <- base_row
        iter_row[names(inipar)] <- model$par
        iter_row["likelihood"] <- model$value
        return(iter_row)
      }

      iter_df <- as.data.frame(do.call(rbind, optim_remaining_pars))
      best_model <- which.min(iter_df$likelihood)
      best_model_AIC <- 2 * iter_df$likelihood[best_model] + 2 * (i + length(basic_pars) + 1)

      print(paste0("iteration: ", iter, ", num parameters included: ", i))
      print(best_model_AIC)
      
      if (best$AIC == 0 | best_model_AIC < best$AIC) # keep the parameter values
        {
          best$AIC <- best_model_AIC
          best$par <- iter_df[best_model, names(all_pars_iter)]
          best$par <- Filter(function(x) !is.na(x), best$par)

          taken <- which(sapply(data_pars_iter, function(x) any(grepl(x, names(best$par)))))

        } else {
          print("No improvement. Exiting loop.")
          break
        }
      
    } # end for i in 1:length(data_pars_iter)
    
  ideal_par_combination <- append(ideal_par_combination, list(best$par))

  write_rds(ideal_par_combination, paste0("./data/", name, "_ideal_par_combination.rds"))

  } #end for iter in 1:nrow(iterations)

  return(ideal_par_combination)

}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------------------- K-fold cross-validation ---------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cross_valid <- function(data, run_function, pars_iter, conditions = NULL) {

  r_squared_values <- c()

  indices <- sample(c(1:5), nrow(data), replace = TRUE)

  for (i in 1:5) {
    # print(i)
    # Define the test and train sets
    test_data <- data[indices == i, ]
    train_data <- data[!indices == i, ]

    # Run the function on the training data and evaluate on the test data
    if (identical(run_function, run_optim)) {
      model_output <- run_function("nls", pars_iter, train_data, conditions, test_data)
    } else {
      model_output <- run_function(train_data, test_data, pars_iter)
    }

    if (is.na(model_output$rsq)) {
      print(paste0(names(pars_iter), ":", model_output$pars))
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
    
    # iter = 1
    # iterations <- iterations_optim
    # print(iter)
    # run_function <- run_optim
    
    
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
      pars_iter <- initial_pars[[iter]] # after definition 
      pars_names <- data_pars_names[[j]]
      
      k <- iterations$basic_par[iter]
      basic_pars_names <- basic_pars_names[[k]]
    } else {
      pars_iter <- data_pars_lm[[j]]
      pars_names <- data_pars_lm_names[[j]]

      basic_pars_iter <- NULL
      basic_pars_names <- NULL
    }

    cross_valid_output <- cross_valid(data, run_function, pars_iter, conditions)

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
  if (export_intermediaries) {
    write.csv(df, paste0("./data/", name, "_results_", model_type, ".csv"), row.names = FALSE)
  }
  
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
    "data_name", "data_pars", "basic_pars", "model_type",
    "rsq", "rsq_sd", "mature_biomass", "age"
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
