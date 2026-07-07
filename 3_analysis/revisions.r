
library(foreach)
library(doParallel)
library(tidyverse)
library(xtable)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_error_propagate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------- One lag per ecoregion of the Amazon ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


data <- import_data("grid_10k_amazon_removed_by_age", biome = 1, n_samples = 30000, asymptote = "nearest_mature", categorical = categorical)

basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(data))[["all"]]

error_prop_results <- error_prop(data, basic_pars, data_pars, conditions, 1000)

results <- data.frame(
    mean_r2 = mean(error_prop_results[[1]]),
    sd_r2 = sd(error_prop_results[[1]]),
    mean_lag = mean(error_prop_results$pars[["lag"]]),
    sd_lag = sd(error_prop_results$pars[["lag"]])
)


write_csv(results, "./0_results/r2_full_amazon_error_prop_by_age.csv")
