# main.R - Main script
# Runs experiments and saves results for different parameter combinations

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(cluster)

# Source other scripts
source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_normalize_cross_validate.r")
source("3_modelling/2_feature_selection_ga.R")

# Set up parallel processing
set.seed(1)
ncore = 4
registerDoParallel(cores = ncore)

biome = 1
n_samples = 10000


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check Importance of parameters included
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# biome = Amazon
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

results <- data.frame()
data <- import_data("grid_10k_amazon_secondary", biome = biome, n_samples = 10000)

# sink("./0_results/amazon_experiments_output.txt")

for (data_pars_name in names(data_pars_options(colnames(data)))) {
    print(data_pars_options(colnames(data))[[data_pars_name]])
    print("------------------------------------------------")

    for (basic_pars_name in names(basic_pars_options)) {
        print(basic_pars_name)

        # Get parameters
        basic_pars <- basic_pars_options[[basic_pars_name]]
        data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

        # Run cross-validation
        cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

        # Return summary
        result <- data.frame(
            basic_pars_name = basic_pars_name,
            data_pars_name = data_pars_name,
            biome = biome,
            mean_r2 = mean(cv_results),
            sd_r2 = sd(cv_results)
        )

        print(result)
        results <- rbind(results, result)
        write.csv(results, file = "./0_results/amazon_experiments.csv", row.names = FALSE)
    }
}

# sink()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save trained model for the best parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# biome = Amazon
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Import data and remove missing values



data <- import_data("grid_10k_amazon_secondary", biome = biome, n_samples = 10000)
head(data)
norm_data <- normalize_independently(data)
# saveRDS(norm_data$train_stats, file = "./0_results/grid_1k_amazon_secondary_train_stats.rds")
norm_data <- norm_data$train_data
head(norm_data)


basic_pars <- basic_pars_options[[basic_pars_name]]
# data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
model <- run_optim(norm_data, init_pars, conditions)
# saveRDS(model, file = paste0("./0_results/amazon_model_", basic_pars_name, ".rds", sep = ""))

 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# compare mean climate and yearly climate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


data <- import_data("grid_10k_amazon_secondary", biome = biome, n_samples = 10000, asymptote = "quarter_biomass")
head(data)
norm_data <- normalize_independently(data)
# saveRDS(norm_data$train_stats, file = "./0_results/grid_1k_amazon_secondary_train_stats.rds")
norm_data <- norm_data$train_data
head(norm_data)

climatic_pars <- c("srad", "temp", "def", "vpd", "pr", "pdsi", "aet")
data_pars_options <- list(
    # mean_climate = colnames(data)[grepl(paste(c(paste0("mean_", climatic_pars)), collapse = "|"), colnames(data))],
    yearly_climate = climatic_pars
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



# use mean climate for the estimated years of growth only





