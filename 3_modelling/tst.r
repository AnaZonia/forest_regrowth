# main.R - Main script
# Runs experiments and saves results for different parameter combinations

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)

# Source other scripts
source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_normalize_cross_validate.r")
source("3_modelling/2_feature_selection.R")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# quick R2 check
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("grid_10k_amazon_secondary_2", biome = 1, n_samples = 10000, asymptote = "nearest_mature")
norm_data <- normalize_independently(data)
norm_data <- norm_data$train_data

basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(norm_data))[["all_mean_climate"]]

init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
model <- run_optim(norm_data, init_pars, conditions)
pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
r2 <- calc_r2(norm_data, pred)
print(r2)
model




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Mean and Yearly Climate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 10000, asymptote = "nearest_mature")
norm_data <- normalize_independently(data)
norm_data <- norm_data$train_data

climatic_pars <- c("srad", "temp", "def", "vpd", "pr", "pdsi", "aet")
data_pars_options <- list(
    mean_climate = colnames(data)[grepl(paste(c(paste0("mean_", climatic_pars)), collapse = "|"), colnames(data))]
    # yearly_climate = climatic_pars
)


for (data_pars_name in names(data_pars_options)) {
    print(data_pars_name)
    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options[[data_pars_name]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
    r2 <- calc_r2(norm_data, pred)
    print(r2)
}

for (data_pars_name in names(data_pars_options)) {
    print(data_pars_name)
    basic_pars <- basic_pars_options[["intercept"]]
    data_pars <- data_pars_options[[data_pars_name]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    pred <- growth_curve(model$par, norm_data)
    r2 <- calc_r2(norm_data, pred)
    print(r2)
}
