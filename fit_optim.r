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
library(tidyverse)
library(fastDummies)
library(foreach)
library(doParallel)

#------------------ SWITCHES ------------------#

run_all <- TRUE
run_one <- TRUE
include_k_par <- FALSE

#------------------ GLobal Variables ------------------#

climatic_vars <- c("prec", "si")

#------------------ FUNCTIONS ------------------#

# - Imports the dataframe
# - Removes unnecessary columns that will not be used in analysis
# - Converts categorical data to dummy variables


# Main Functions

import_data <- function(path) {
  data <- read_csv(path) %>%
    select(-c(.geo, latitude, longitude)) %>%
    select(-starts_with("system"))

  # Convert specified variables to factors
  categorical <- c("ecoreg", "soil", "last_LU")
  # Convert categorical variables to factors
  data <- data %>%
    mutate(across(all_of(categorical), as.factor))
  
  # data <- data %>%
  #   group_by(ecoreg) %>%
  #   mutate(
  #     mean_per_ecoregion = mean(mature_biomass, na.rm = TRUE),
  #     sd_per_ecoregion = sd(mature_biomass, na.rm = TRUE)
  #   ) %>%
  #   ungroup()

  # data <- data %>%
  #   filter(
  #     mature_biomass > agbd,
  #     mature_biomass < (mean_per_ecoregion + sd_per_ecoregion),
  #     mature_biomass > (mean_per_ecoregion - sd_per_ecoregion)
  #   )

  # Create dummy variables
  data <- dummy_cols(data, select_columns = categorical, remove_selected_columns = TRUE)
  data <- data %>% filter(!is.na(mature_biomass))

  data
}

import_climatic_data <- function(path, normalize) {
  data <- import_data(path)

  means <- sapply(climatic_vars, function(var) rowMeans(data[, grep(var, names(data))], na.rm = TRUE))
  colnames(means) <- paste0("mean_", climatic_vars)
  data <- cbind(data, means)

  df_climatic_hist <- tibble()
  for (age in 1:max(data$age)) {
    age_data <- data %>% filter(age == .env$age)
    years <- seq(2019, 2019 - age + 1, by = -1)
    # Identify all columns including the variables in climatic_vars
    clim_columns <- expand.grid(climatic_vars, years) %>%
      unite(col = "col", sep = "_") %>%
      pull(col)

    # subsect the dataframe to only include the climatic columns of the desired years
    all_clim_columns <- names(data)[str_detect(names(data), paste(climatic_vars, "_", collapse = "|"))]

    # turn all values in the columns of years not included to 0
    clim_columns_not_included <- setdiff(all_clim_columns, clim_columns)
    age_data[clim_columns_not_included] <- 0

    df_climatic_hist <- bind_rows(df_climatic_hist, age_data)
  }

  if (normalize) {
    df_climatic_hist <- df_climatic_hist %>%
      mutate(across(
        where(is.numeric) &
          !matches("soil|biome|ecoreg|last_LU|protec|indig|agbd|mature_biomass|mean_per_ecoregion|sd_per_ecoregion"),
        ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
      )) %>%
      select(where(~ sum(is.na(.)) < nrow(df_climatic_hist)))
  }
  df_climatic_hist
}


# Intakes pars <- a vector with the initial parameters to be included
# data <- the dataframe with the predictors
# pars_chosen <- the list of parameters to be added into the shape term

# An example would look like:
#   growth_curve(c(B0 = 40, A = 80, theta = 5, age = 2), data_10[[1]], c('age'))

# define the climatic parameters - the ones that change yearly
growth_curve <- function(pars, data, pars_chosen) {
  k <- rep(0, nrow(data))

  if (include_k_par) {
    k <- k + pars[["k"]]
  }

  non_clim_vars <- setdiff(pars_chosen, climatic_vars)

  tmp_age <- 1
  tmp_age1 <- data[[age]]
  if ("age" %in% pars_chosen)
  {
    tmp_age <- data[["age"]]
    tmp_age1 <- 1
  }

  #if (any(climatic_vars %in% pars_chosen)) {
  k <- k + rowSums(sapply(non_clim_vars, function(var) pars[[var]] * data[[var]] * tmp_age))

  for (clim_var in intersect(climatic_vars, pars_chosen)) {
    years <- seq(2019, 1985, by = -1)
    clim_columns <- paste0(clim_var, "_", years)
    k <- k + rowSums(sapply(clim_columns, function(col) pars[[clim_var]] * data[[col]/tmp_age1]))
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
  # "data/5y_LULC.csv",
  "data/10y_LULC_1km.csv" # "data/15y_LULC_1km.csv",
  # "data/all_LULC.csv"
)



dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

# adding names can make indexing complicated, so they are being referenced as a vector
# names_dataframes <- c("data_5", "data_10", "data_15", "data_all")
# dataframes <- lapply(dataframes, function(df) {
#   df %>% filter(!is.na(mature_biomass))
# })


if (run_all || run_one) {
  # Define conditions
  conditions <- list(
    'pars["theta"] > 10',
    'pars["theta"] < 0',
    'pars["B0"] < 0'
  )

  # intercept, shape term, growth rate
  pars_basic <- c(B0 = 40, theta = 5)
  if (include_k_par) {
    pars_basic <- c(pars_basic, k = 0.01)
  }

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
    c("indig", "protec", "cwd", "mean_prec", "mean_si", names(data)[str_detect(names(data), "LU")]),
    colnames_filtered,
    c(colnames_filtered, climatic_vars)
  )

  # adding names can make indexing complicated, so they are being referenced as a vector
  names_configurations <- c("age", "fires", "age_fires", "all_cat", "all_non_hist", "all_hist")
}


if (run_one) {
  data <- dataframes[[1]]
  data <- data %>%
    filter(
      mature_biomass > agbd
    )
  pars_chosen <- configurations[[4]]
  val <- run_optimization(
    "nls", pars_basic, data, pars_chosen,
    if ("age" %in% pars_chosen) {
      c(conditions, list('pars["age"] < 0', 'pars["age"] > 5'))
    } else {
      conditions
    }
  )
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

data1 <- dataframes[[1]] %>% filter(biome == 1)
plot(agbd ~ mature_biomass, data = data1, col = "blue")
abline(0, 1, col = "red")

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
      pars[["B0"]] + (pars[["A"]] - pars[["B0"]]) * (1 - exp(-pars[["k"]] * data$age))^pars[["theta"]]
    }
  } else {
    function(pars, data) {
      pars[["B0"]] + (data[["mature_biomass"]] - pars[["B0"]]) * (1 - exp(-pars[["k"]] * data$age))^pars[["theta"]]
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
  # print(names_dataframes[i])

  # Without asymptote (using mature_biomass)
  result1 <- run_growth_model(dataframes[[i]], c(B0 = 40, theta = 5, k = 0.1))
  print(paste("R-squared (fixed asymptote, fit growth rate):", result1$r_squared))

  # With asymptote
  result2 <- run_growth_model(dataframes[[i]], c(B0 = 40, theta = 5, k = 0.1, A = 100))
  print(paste("R-squared (fit asymptote, rate fit from age):", result2$r_squared))
}



library(mgcv)
library(randomForest)

data <- dataframes[[1]]
# data <- read.csv("data/10y_LULC_1km.csv")
# data <- data %>% filter(!is.na(mature_biomass))
# nrow(data)
# mean(data$agbd)
# mean(data$mature_biomass)

# lulc_columns <- names(data)[str_detect(names(data), "lulc")]
# lulc_unique_values <- lapply(lulc_columns, function(col) unique(data[[col]]))
# names(lulc_unique_values) <- lulc_columns
# # Print the results
# for (col in lulc_columns) {
#   cat(col, ":\n")
#   print(lulc_unique_values[[col]])
#   cat("\n")
# }
# Split data into training and testing sets

set.seed(123)
train_indices <- sample(1:nrow(data), size = floor(0.7 * nrow(data)))
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]


lulc_sum_columns <- names(data)[str_detect(names(data), "lulc")]
lulc_sum_columns <- lulc_sum_columns[!grepl("lulc_sum_(48|47|40)$", lulc_sum_columns)]

pars_categ <- c("indig", "protec", names(data)[str_detect(names(data), "LU")])
# pars_categ <- c()
pars_smooth <- c("mean_si", "mean_prec", "cwd", "mature_biomass"
# , lulc_sum_columns # "lulc_sum_15"
)

# Fit a GAM model
formula <- as.formula(paste(
  "agbd ~",
  paste(sapply(pars_smooth, function(x) paste0("s(", x, ")")), collapse = " + "),
  "+",
  paste(pars_categ, collapse = " + ")
))

gam_model <- gam(formula, data = train_data)
summary(gam_model)

# Predict using the GAM model
pred_gam <- predict(gam_model, newdata = test_data)
cor(test_data$agbd, pred_gam)^2

# Plot predictions vs observed
plot(pred_gam, data$agbd)
abline(0, 1, col = "red")

rf_pars <- c(pars_smooth, pars_categ, "age")
rf_pars

rf_formula <- as.formula(paste("agbd ~", paste(rf_pars, collapse = " + ")))
summary(lm(rf_formula, data = data))
summary(lm(agbd ~ mature_biomass, data = data))

train_data <- na.omit(train_data)

# Fit a Random Forest model on the training data
rf_model <- randomForest(rf_formula, data = train_data, ntree = 500, mtry = sqrt(ncol(train_data) - 1), importance = TRUE, keep.forest = TRUE, oob.score = TRUE)

# Print model summary
print(rf_model)

# Predict using the Random Forest model on the test data
pred_rf_test <- predict(rf_model, newdata = test_data)

calc_r_squared(test_data$agbd, pred_rf_test)

# Plot predictions vs observed for the test data
plot(pred_rf_test, test_data$agbd)
abline(0, 1, col = "red")



data <- read.csv("./data/mat_for_biomass_and_distance.csv")
data <- data[, 2:3]

# why in here it is showing only a few points so far away? makes no sense.
data <- subset(data, distance < 5000)
plot(data$distance, data$mature_biomass)
