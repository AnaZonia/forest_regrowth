library(mgcv)
library(randomForest)
library(tidyverse)

#------------------ FUNCTIONS ------------------#


# Intakes pars <- a vector with the initial parameters to be included
# data <- the dataframe with the predictors
# pars <- the list of parameters to be added into the shape term

growth_curve <- function(pars, data) {

  non_clim_pars <- setdiff(names(pars), c(non_data_pars, climatic_pars))

  implicit_age <- if (!"age" %in% names(pars)) data[["age"]] else rep(1, nrow(data))
  
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

run_optimization <- function(fun, pars_basic, pars_chosen, train_data, test_data, conditions) {
  # pars_chosen <- data_pars[[2]]
  # pars_basic <- basic_pars[[4]]
  # train_data <- train_dataframes[[1]]
  # test_data <- test_dataframes[[1]]
  # fun = "nls"
  # Run optimization
  if (fun == "nll") {
    pars <- c(pars, setNames(0.1, "sd"))
  }

  current_pars <- c(
    theta = 5,
    setNames(
      rep(0.1, length(c(pars_basic, pars_chosen[1]))),
      c(pars_basic, pars_chosen[1])
    )
  )
  current_pars
  if ("B0" %in% pars_basic) {
    current_pars["B0"] <- 40
  }
  
  # Function to perform optimization
  optimize <- function(pars) {
    optim(pars,
      likelihood,
      fun = fun,
      data = train_data,
      conditions = conditions
    )
  }

  o <- optimize(current_pars)

  if (length(pars_chosen) > 1) {
    for (i in 2:length(pars_chosen)) {
      current_pars <- c(o$par, setNames(0.1, pars_chosen[i]))
      o <- optimize(current_pars)
    }
  }

  pred <- growth_curve(o$par, test_data)
  rsq <- cor(test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))
  return(results <- list(
    model = o,
    rsq = rsq
  ))
}



run_gam_lm <- function(train_data, test_data, pars_chosen, model) {

  # Fit a GAM model
  if (model == "gam") {
    # Filter data_pars to exclude items containing climatic_pars
    continuous <- c(pars_chosen[!pars_chosen %in% pars_categ])
    categorical <- pars_chosen[pars_chosen %in% pars_categ]

    formula <- as.formula(paste(
      "agbd ~",
      paste(sapply(continuous, function(x) paste0("s(", x, ")")), collapse = " + "),
      "+",
      paste(categorical, collapse = " + ")
    ))

    model <- gam(formula, data = train_data)
  } else { # Fit a Linear Model

    rf_lm_formula <- as.formula(paste("agbd ~", paste(pars_chosen, collapse = " + ")))

    if (model == "lm") {
      model <- lm(rf_lm_formula, data = train_data)
      # Check for rank deficiency
      aliased_vars <- summary(model)$aliased

      if (any(aliased_vars)) {
        problematic_vars <- names(aliased_vars)[aliased_vars]
        print(paste("Removing rank-deficient variables:", paste(problematic_vars, collapse = ", ")))

        # Remove problematic variables from the formula
        pars_chosen <- pars_chosen[!pars_chosen %in% problematic_vars]
        rf_lm_formula <- as.formula(paste("agbd ~", paste(pars_chosen, collapse = " + ")))
        model <- lm(rf_lm_formula, data = train_data)
      }
    }
  }

  pred <- predict(model, newdata = test_data)
  rsq <- cor(test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))

  return(results <- list(
    model = model,
    rsq = rsq
  ))
}


run_rf <- function(train_data, test_data, pars_chosen) {
  rf_lm_formula <- as.formula(paste("agbd ~", paste(pars_chosen, collapse = " + ")))

  sampled_train_data <- train_data[sample(nrow(train_data), 20000), ]

  model <- randomForest(rf_lm_formula,
    data = sampled_train_data,
    ntree = 100, mtry = 2, importance = TRUE,
    keep.forest = TRUE, oob.score = TRUE, do.trace = 10, parallel = TRUE
  )

  pred <- predict(model, newdata = test_data)
  rsq <- cor(test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))
  return(results <- list(
    output = importance(rf_iter$model),
    rsq = rsq
  ))
}



# Define helper functions
process_row <- function(output, model_name, model_type, rsq) {
  row <- as.data.frame(t(output))
  missing_cols <- setdiff(unique(c(unlist(c(data_pars, non_data_pars)), "mature_biomass", "age")), names(row))
  row[missing_cols] <- NA
  row <- row[, unique(c(unlist(c(data_pars, non_data_pars)), "mature_biomass", "age"))]
  row$model_name <- model_name
  row$model_type <- model_type
  row$rsq <- rsq
  row %>% select(model_name, model_type, rsq, "mature_biomass", "age", all_of(non_data_pars), everything())
}
