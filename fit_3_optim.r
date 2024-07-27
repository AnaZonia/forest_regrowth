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
set.seed(1)

#------------------ SWITCHES ------------------#

run_all <- TRUE
run_one <- FALSE

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

run_optimization <- function(fun, pars_basic, data_pars, train_data, test_data, conditions) {

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
      data = train_data,
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

  pred <- growth_curve(o$par, test_data)
  rsq <- cor(test_data$agbd, pred)^2
  print(paste("R-squared:", rsq))
  return(results <- list(
    model = o,
    rsq = rsq
  ))
}

run_rf_gam_lm <- function(train_data, test_data, categorical, continuous, model) {

  # Fit a GAM model
  if (model == "gam") {
    formula <- as.formula(paste(
      "agbd ~",
      paste(sapply(continuous, function(x) paste0("s(", x, ")")), collapse = " + "),
      "+",
      paste(categorical, collapse = " + ")
    ))
    model <- gam(formula, data = train_data)
  } else {
    rf_lm_pars <- c(continuous, categorical)
    rf_lm_formula <- as.formula(paste("agbd ~", paste(rf_lm_pars, collapse = " + ")))
    if (model == "lm") {
      model <- lm(rf_lm_formula, data = train_data)
      # Check for rank deficiency
      aliased_vars <- summary(lm_iter$model)$aliased

      if (any(aliased_vars)) {
        problematic_vars <- names(aliased_vars)[aliased_vars]
        print(paste("Removing rank-deficient variables:", paste(problematic_vars, collapse = ", ")))

        # Remove problematic variables from the formula
        rf_lm_pars <- rf_lm_pars[!rf_lm_pars %in% problematic_vars]
        rf_lm_formula <- as.formula(paste("agbd ~", paste(rf_lm_pars, collapse = " + ")))
        model <- lm(rf_lm_formula, data = train_data)
      }

    } else {
      model <- randomForest(rf_lm_formula,
        data = train_data,
        ntree = 500, mtry = sqrt(ncol(train_data) - 1), importance = TRUE,
        keep.forest = TRUE, oob.score = TRUE
      )
    }
  }

  pred <- predict(model, newdata = test_data)
  rsq <- cor(test_data$agbd, pred)^2
  print(rsq)

  return(results <- list(
    model = model,
    pred = pred,
    test_agbd = test_data$agbd,
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

#------------------ Global Variables ------------------#

climatic_pars <- c("prec", "si")
# parameters not fit with data - just for the functional form
non_data_pars <- c("k0", "B0_exp", "B0", "theta")

pars_categ <- c("indig", "protec", names(data)[str_detect(names(data), "LU|ecoreg|soil")])

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

# Function to create both samples
sample_data <- function(df, size_train, size_test) {
  train_dataframes <- sample_n(df, size = size_train)
  remaining <- df[!(row.names(df) %in% row.names(train_dataframes)), ]
  test_dataframes <- sample_n(remaining, size = size_test)
  list(train_dataframes, test_dataframes)
}

# Apply the function to each dataframe
samples <- lapply(dataframes, sample_data, size_train = 60000, size_test = 40000)

# Extract the first and second element of samples
train_dataframes <- lapply(samples, `[[`, 1)
test_dataframes <- lapply(samples, `[[`, 2)

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
    c("cwd"),
    c("cwd", "mean_prec", "mean_si"),
    c("cwd", climatic_pars),
    colnames_filtered[!grepl("ecoreg|soil", colnames_filtered)],
    colnames_filtered,
    c(colnames_filtered, "cwd", "mean_prec", "mean_si"),
    c(colnames_filtered, climatic_pars)
  )

  basic_pars <- list(
    c("age", "B0"),
    c("age", "k0", "B0"),
    c("k0", "B0"), # age is implicit
    c("k0", "B0_exp") # age is implicit
  )
}

if (run_one) {
  pars_chosen <- data_pars[[1]]
  pars_basic <- basic_pars[[1]]
  train_data <- train_dataframes[[1]]
  test_data <- test_dataframes[[1]]

  conditions_iter <- conditions
  if ("age" %in% pars_chosen) {
    conditions_iter <- c(conditions_iter, list('pars["age"] < 0', 'pars["age"] > 5'))
  } else if ("B0" %in% pars_chosen) {
    conditions_iter <- c(conditions_iter, list('pars["B0"] < 0'))
  }
  
  o_iter <- run_optimization("nls", pars_basic, pars_chosen, train_data, test_data, conditions_iter)

}

if (run_all) {
  start_time <- Sys.time()
  print(start_time)

  iterations_optim <- expand.grid(
    dataframe = seq_along(train_dataframes),
    data_par = seq_along(data_pars),
    basic_par = seq_along(basic_pars)
  )

  registerDoParallel(cores = 15)

  results_optim <- foreach(iter = 1:nrow(iterations_optim),  .combine = 'bind_rows', .packages = c('dplyr')) %dopar% {
    print(iter)
    i <- iterations_optim$dataframe[iter]
    j <- iterations_optim$data_par[iter]
    k <- iterations_optim$basic_par[iter]

    train_data <- train_dataframes[[i]]
    test_data <- test_dataframes[[i]]
    pars_chosen <- data_pars[[j]]
    pars_basic <- basic_pars[[k]]
    
    conditions_iter <- conditions
    if ("age" %in% pars_chosen) {
      conditions_iter <- c(conditions_iter, list('pars["age"] < 0', 'pars["age"] > 5'))
    } else if ("B0" %in% pars_chosen) {
      conditions_iter <- c(conditions_iter, list('pars["B0"] < 0'))
    }

    o_iter <- run_optimization("nls", pars_basic, pars_chosen, train_data, test_data, conditions_iter)

    optim_output <- o_iter$model$par

    # Reorder columns
    optim_row <- process_row(optim_output, intervals[[i]], "optim", o_iter$rsq)

    print(optim_row)
    optim_row
  }


  data_pars_lm <- append(
    lapply(data_pars, function(x) c(x, "mature_biomass")),
    list(c("age"), c("mature_biomass"), c("age", "mature_biomass"))
  )
  
  iterations_lm <- expand.grid(
    dataframe = seq_along(dataframes),
    data_par = which(!sapply(data_pars_lm, function(par_set) any(climatic_pars %in% par_set)))
  )

  results_lm <- foreach(iter = 1:nrow(iterations_lm), .combine = "bind_rows", .packages = c("dplyr")) %dopar% {

    i <- iterations_lm$dataframe[iter]
    j <- iterations_lm$data_par[iter]

    train_data <- train_dataframes[[i]]
    test_data <- test_dataframes[[i]]
    pars_chosen <- data_pars_lm[[j]]

    # Filter data_pars to exclude items containing climatic_pars
    continuous <- c(pars_chosen[!pars_chosen %in% pars_categ])
    categorical <- pars_chosen[pars_chosen %in% pars_categ]

    lm_iter <- run_rf_gam_lm(train_data, test_data, categorical, continuous, "lm")

    lm_output <- summary(lm_iter$model)$coefficients[-1, 1, drop = FALSE] # -1 to remove (Intercept)

    lm_row <- process_row(lm_output, intervals[[i]], "lm", lm_iter$rsq)

    print(lm_row)
    lm_row
  }

  # Combine all results into a single dataframe
  df <- as.data.frame(rbind(results_optim, results_lm))

  # Write the dataframe to a CSV file
  write.csv(df, "./data/fit_results_test_train.csv", row.names = FALSE)
  print(paste("written! Time took: ", as.numeric(difftime(Sys.time(), start_time, units = "hours")), " hours"))

}

