
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
ncore <- 4
registerDoParallel(cores = ncore)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save trained model for the best parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000)
norm_data <- normalize_independently(data)

norm_data <- norm_data$train_data

for (basic_pars_name in names(basic_pars_options)) {
    basic_pars <- basic_pars_options[[basic_pars_name]]

    if (basic_pars_name == "intercept") {
        # For intercept, we force the data to intercept through zero
        mean_biomass_at_zero_age <- median(norm_data$biomass[norm_data$age == 1], na.rm = TRUE)
        norm_data$biomass <- norm_data$biomass - mean_biomass_at_zero_age
    }

    data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    saveRDS(model, file = paste0("./0_results/amazon_model_", basic_pars_name, ".rds", sep = ""))
    
    pred_vs_obs <- data.frame(
        age = norm_data$age,
        pred = growth_curve(model$par, norm_data)
    )

    print(paste0("R2 for ", basic_pars_name, ": ", calc_r2(norm_data, pred_vs_obs$pred)))

    # get the biomass predictions for the ages 1 to 105 (in 35 year steps)
    for (i in c(35, 70)) {
        norm_data_future <- norm_data
        norm_data_future$age <- norm_data_future$age + i
        pred_future <- growth_curve(model$par, norm_data_future)
        df_future <- data.frame(
            age = norm_data_future$age,
            pred = pred_future
        )
        pred_vs_obs <- rbind(pred_vs_obs, df_future)
    }
    
    if (basic_pars_name == "intercept") {
        pred_vs_obs$pred <- round(pred_vs_obs$pred - model$par["B0"])
    } else {
        # add column obs with the age correspondent to that in norm_data
        pred_vs_obs <- cbind(pred_vs_obs, data.frame(
            obs = norm_data$biomass,
            obs_age = round(norm_data$age + model$par["lag"])
        ))
    }

    write.csv(pred_vs_obs, file = paste0("./0_results/pred_vs_obs_amazon_", basic_pars_name, ".csv"), row.names = FALSE)

}

