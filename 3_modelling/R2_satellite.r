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
ncore = 4
registerDoParallel(cores = ncore)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Compare parameter combinations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

results <- data.frame()
data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 10000)


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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Biomes (Atlantic, Amazon)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Asymptotes ("nearest_mature", "ecoreg_biomass", "quarter_biomass", "full_amazon")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Land Use (non_aggregated_all, aggregated_all, non_aggregated_5yr, non_aggregated_10yr)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




