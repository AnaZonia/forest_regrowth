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
library(mgcv)
library(randomForest)

source("fit_1_import_data.r")
set.seed(1)

#------------------ SWITCHES ------------------#

run_all <- TRUE

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

################### Running model ###################

# pars_chosen <- data_pars[[1]]
# pars_basic <- basic_pars[[4]]
# train_data <- train_dataframes[[1]]
# test_data <- test_dataframes[[1]]

# conditions_iter <- conditions
# if ("age" %in% pars_chosen) {
#   conditions_iter <- c(conditions_iter, list('pars["age"] < 0', 'pars["age"] > 5'))
# } else if ("B0" %in% pars_chosen) {
#   conditions_iter <- c(conditions_iter, list('pars["B0"] < 0'))
# }

# o_iter <- run_optimization("nls", pars_basic, pars_chosen, train_data, test_data, conditions_iter)


if (run_all) {
  start_time <- Sys.time()
  print(start_time)

  iterations_optim <- expand.grid(
    dataframe = seq_along(train_dataframes),
    data_par = seq_along(data_pars),
    basic_par = seq_along(basic_pars)
  )

  registerDoParallel(cores = 25)

  results_optim <- foreach(iter = 1:nrow(iterations_optim), .combine = 'bind_rows', 
  .packages = c('dplyr')) %dopar% {
    
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
    print(paste("Time so far: ", as.numeric(difftime(Sys.time(), start_time, units = "mins")), " minutes"))
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
    print(iter)
    i <- iterations_lm$dataframe[iter]
    j <- iterations_lm$data_par[iter]

    train_data <- train_dataframes[[i]]
    test_data <- test_dataframes[[i]]
    pars_chosen <- data_pars_lm[[j]]

    lm_iter <- run_gam_lm(train_data, test_data, pars_chosen, "lm")

    lm_output <- summary(lm_iter$model)$coefficients[-1, 1, drop = FALSE] # -1 to remove (Intercept)

    lm_row <- process_row(lm_output, intervals[[i]], "lm", lm_iter$rsq)

    print(lm_row)
    print(paste("Time so far: ", as.numeric(difftime(Sys.time(), start_time, units = "mins")), " minutes"))
    lm_row
  }


  results_rf <- foreach(
    iter = 1:nrow(iterations_lm), .combine = "bind_rows",
    .packages = c("dplyr", "randomForest")
  ) %dopar% {

    print(iter)
    i <- iterations_lm$dataframe[iter]
    j <- iterations_lm$data_par[iter]

    train_data <- train_dataframes[[i]]
    test_data <- test_dataframes[[i]]
    pars_chosen <- data_pars_lm[[j]]

    rf_iter <- run_rf(train_data, test_data, pars_chosen)

    # writing only %IncMSE into the dataframe
    rf_row <- process_row(rf_iter$output[, 1], intervals[[i]], "rf", rf_iter$rsq)
    
    print(paste("Time so far: ", as.numeric(difftime(Sys.time(), start_time, units = "mins")), " minutes"))
    print(rf_row)
    rf_row
  }
  
  
  print(paste("written! Time for the whole operation: ", as.numeric(difftime(Sys.time(), start_time, units = "mins")), " minutes"))


  # Combine all results into a single dataframe
  df <- as.data.frame(rbind(results_optim, results_lm, results_rf))

  # Write the dataframe to a CSV file
  write.csv(df, "./data/fit_results_test_train_rf.csv", row.names = FALSE)
  print(paste("written! Time took: ", as.numeric(difftime(Sys.time(), start_time, units = "hours")), " hours"))

}

