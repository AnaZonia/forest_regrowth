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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("./6_old_scripts/fit_1_import_data.r")

# data <- import_data("./0_data/data_mapbiomas_4.csv")
data <- import_data("./0_data/non_aggregated.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define different sets of data parameters for modeling
data_pars <- list(
  c("nitro", "sur_cover"),
  c("nitro", "sur_cover", "lulc_sum_15", "num_fires_before_regrowth", "distance", "fallow"),
  names(data)[!names(data) %in% c("agbd", "nearest_mature_biomass", "age")]
)

data_pars_names <- c(
  "nitro", "nitrosur", "most", "all"
)


data <- read_csv("./0_data/non_aggregated.csv", show_col_types = FALSE)
data <- subset(data, biome == 4)

data <- data[, !names(data) %in% c(
  "biome", "ecoreg", "topography",
  "cwd", "mean_aet", "lulc_sum_20", "lulc_sum_40", "lulc_sum_48"
)]

categorical <- c("last_LU", "indig", "protec")
data <- data %>%
  mutate(across(all_of(categorical), as.factor))
data <- data %>%
  filter(!last_LU %in% c(20, 40, 48))

# vars_to_keep <- c(
#   names(data)[str_detect(names(data), "lulc_sum")],
#   "last_LU",
#   "fallow",
#   "num_fires_before_regrowth",
#   "indig",
#   "protec"
# )
vars_to_keep <- names(data)[!names(data) %in% c("agbd", "nearest_mature_biomass", "age")]

data <- data %>%
  mutate(across(
    where(is.numeric) &
      matches(paste(vars_to_keep, collapse = "|")),
    ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
  )) %>%
  select(where(~ sum(is.na(.)) < nrow(data)))


formula <- as.formula(paste("agbd ~", paste(vars_to_keep, collapse = " + ")))
summary(lm(formula, data))


vars_to_keep


# Define basic parameter sets for modeling
basic_pars <- list(
  c("age", "B0"),
  c("age", "k0", "B0"),
  c("k0", "B0"), # age is implicit
  c("k0", "B0_exp") # age is implicit
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
  initial_pars <- find_combination_pars(iterations_optim, "nls")
} else {
  initial_pars <- readRDS(paste0("./0_results/", name, "_ideal_par_combination.rds"))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Run Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# test_rows <- which(iterations_optim$interval == 1 & iterations_optim$data_par == 5 &
#   iterations_optim$basic_par == 1)

start_time <- Sys.time()
print(start_time)

results_optim <- run_foreach(iterations_optim, "optim", run_optim, conditions)

results_all <- results_optim %>%
  arrange(data_pars, basic_pars, data_name) # sort rows by parameter

write.csv(results_all, paste0("./data/", name, "_results_all.csv"), row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

print(paste(
  region, "ALL DONE! Time for the whole thing using ", ncores, "cores: ",
  as.numeric(difftime(Sys.time(), start_time, units = "hours")), "hours"
))
