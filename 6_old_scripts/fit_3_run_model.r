library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)

# Source external R scripts for data import and function definitions
source("./6_old_scripts/fit_1_import_data.r")
source("./6_old_scripts/fit_2_functions.r")

set.seed(1)
ncores <- 30
registerDoParallel(cores = ncores)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------------- Switches ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

test_switch <- FALSE # Set to TRUE for test_rows or FALSE to run all rows
# include test_rows rows to test
test_rows <- 1

run_initial_fit_by_order <- TRUE
export_intermediaries <- TRUE

region <- "atla"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# number of rows to be included in analysis
n_samples <- 10000

# name for export
name <- region

# Define conditions for parameter constraints
conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0')
non_data_pars <- c("k0", "B0_exp", "B0", "theta")

climatic_pars <- c("prec", "si")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("./0_data/non_aggregated_all.csv")
# data <- import_data("./data/non_aggregated_all.csv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define different sets of data parameters for modeling
data_pars <- list(
    c("sur_cover"),
    # c("sur_cover", "lulc_sum_15", "num_fires_before_regrowth", "distance", "fallow"),
    # c("sur_cover", "num_fires_before_regrowth", "distance", "fallow", "last_LU_9", "last_LU_15", "last_LU_20", "last_LU_21", "last_LU_39", "last_LU_41", "last_LU_46", "last_LU_48", "lulc_sum_9", "lulc_sum_15", "lulc_sum_20", "lulc_sum_21", "lulc_sum_39", "lulc_sum_40", "lulc_sum_41", "lulc_sum_46", "lulc_sum_48", "protec", "indig")
    names(data)[!names(data) %in% c("agbd", "nearest_mature_biomass", "age")]
)

data_pars_names <- c(
    "sur", "all"
)

# Define basic parameter sets for modeling
basic_pars <- list(
    c("age", "B0", "theta"),
    c("age", "k0", "B0", "theta"),
    c("k0", "B0", "theta") # age is implicit
)

basic_pars_names <- as.list(sapply(basic_pars, function(par_set) paste(par_set, collapse = "_")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create grids of different combinations of model inputs

# Optim / growth curve
iterations_optim <- expand.grid(
    data_par = seq_along(data_pars),
    basic_par = seq_along(basic_pars)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------------------- Find ideal parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if (run_initial_fit_by_order) {
    # keep combination of parameters for initial fit
    initial_pars <- find_combination_pars(iterations_optim, conditions)
} else {
    initial_pars <- readRDS(paste0("./0_data/", name, "_ideal_par_combination.rds"))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Run Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# test_rows <- which(iterations_optim$interval == 1 & iterations_optim$data_par == 5 &
#   iterations_optim$basic_par == 1)
# source("./6_old_scripts/fit_2_functions.r")
# initial_pars <- readRDS(paste0("./0_data/atla_ideal_par_combination.rds"))
initial_pars <- find_combination_pars(iterations_optim, conditions)

start_time <- Sys.time()
print(start_time)

source("./6_old_scripts/fit_2_functions.r")
results_optim <- run_foreach(iterations_optim, conditions)
results_optim$rsq


results_all <- results %>%
  arrange(data_pars, basic_pars, data_name) # sort rows by parameter

write.csv(results_all, paste0("./0_data/", name, "_results_all.csv"), row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Define a function to run foreach for different models
run_foreach <- function(iterations, model_type, run_function, conditions = NULL) {
  results <- foreach(
    iter = if (test_switch) test_rows else 1:nrow(iterations),
    .combine = "bind_rows", .packages = c("dplyr", "randomForest")
  ) %dopar% {
    i <- iterations$interval[iter]
    j <- iterations$data_par[iter]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

    pars_iter <- initial_pars[[iter]] # after definition
    pars_names <- data_pars_names[[j]]

    k <- iterations$basic_par[iter]
    basic_pars_names <- basic_pars_names[[k]]

    cross_valid_output <- cross_valid(data, run_function, pars_iter, conditions)

    # Organize result
    row <- process_row(cross_valid_output, model_type, pars_names,
      basic_pars_names = basic_pars_names
    )

    print(paste("Time so far: ", as.numeric(difftime(
      Sys.time(), start_time,
      units = "mins"
    )), " minutes"))
    print(row)

    row
  }

  # Calculate and print the total time taken
  print(paste(
    model_type, "finished! Time for the whole operation: ",
    as.numeric(difftime(Sys.time(), start_time, units = "hours")), " hours"
  ))


  write.csv(df, paste0("./0_results/", name, "_results_optim.csv"), row.names = FALSE)

  return(df)
}
