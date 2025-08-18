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
source("3_modelling/2_cross_validate.r")
source("3_modelling/2_feature_selection.r")

# Set up parallel processing
set.seed(1)
ncore = 4
registerDoParallel(cores = ncore)

data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000, asymptote = "nearest_mature")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Asymptotes ("nearest_mature", "ecoreg_biomass", "quarter_biomass", "full_amazon")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
results <- data.frame()

for (asymptote in c("nearest_mature", "ecoreg_biomass", "quarter_biomass", "full_amazon")) {
    data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000, asymptote = asymptote)

    data_pars_name <- "age_only"
    basic_pars_name <- "lag"

    # Get parameters
    basic_pars <- basic_pars_options[[basic_pars_name]]
    data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

    # Run cross-validation
    cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

    # Return summary
    result <- data.frame(
        basic_pars_name = basic_pars_name,
        asymptote = asymptote,
        mean_r2 = mean(cv_results),
        sd_r2 = sd(cv_results)
    )

    print(result)
    results <- rbind(results, result)
    write.csv(results, file = "./0_results/R2_asymptotes.csv", row.names = FALSE)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Land Use (non_aggregated_all, aggregated_all, non_aggregated_5yr, non_aggregated_10yr)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

results <- data.frame()

land_use_list <- list.files(paste0("./0_data"), pattern = "land_use", full.names = FALSE)
land_use_list

for (biome in c(1, 4)) {
    print(biome)

    data <- import_data("grid_10k_secondary_non_aggregated_5yr", biome_num = biome, n_samples = 10000)
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
            write.csv(results, file = "./0_results/R2_satellite_5yr.csv", row.names = FALSE)
        }
    }
}






asymptote <- "nearest_mature"
data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000, asymptote = asymptote)

data_pars_name <- "all_mean_climate"
basic_pars_name <- "lag"

# Get parameters
basic_pars <- basic_pars_options[[basic_pars_name]]
data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

data_pars <- data_pars[!data_pars %in% c("mean_def")]

# Run cross-validation
cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

# Normalize training and test sets independently, but using training data's min/max for both
train_data <- normalize_independently(data)$train_data

# Function to perform direct optimization
pars_init <- find_combination_pars(basic_pars, data_pars, train_data)

# Run the model function on the training set and evaluate on the test set
model <- run_optim(train_data, pars_init, conditions)

pred_cv <- growth_curve(model$par, train_data, lag = model$par["lag"])

# save the predicted values of each iteration of the cross validation.
r2 <- calc_r2(train_data, pred_cv)
r2
