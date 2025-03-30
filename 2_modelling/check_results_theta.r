
# Load necessary libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

# get data for plotting
source("2_modelling/1_modelling.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/1_parameters.r")

biome <- 1
n_samples <- 10000
set.seed(1)
registerDoParallel(cores = 4)

conditions <- list('pars["k0"] < 0')

thetas <- c(0.5, 1, 2, NA)
raw_age <- c(TRUE, FALSE)
fit_method <- c("B0", "lag")

df <- data.frame("raw_age" = NA, "fit_method" = NA, "theta" = NA, "k0" = NA, "B0" = NA,  "lag" = NA, "R2" = NA)

for (theta in thetas) {
    for (raw_age in raw_age) {
        for (fit_method in fit_method) {
            
            for (i in 1:5) {
            
                data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples) %>%
                    select(-mean_pdsi)

                # Fit the model on the full data
                norm <- normalize_independently(data, raw_age)
                norm_data <- norm$train_data
                norm_stats <- norm$train_stats

                head(norm_data[, c("age", "num_fires", "sur_cover", "nearest_biomass", "biomass")])

                basic_pars <- c("k0", fit_method)

                if (is.na(theta)) {
                    basic_pars <- c(basic_pars, "theta")
                }

                pars_init <- find_combination_pars(basic_pars, "data_pars" = c("num_fires", "sur_cover"), norm_data)

                final_model <- run_optim(norm_data, pars_init, conditions)

                if (!is.na(final_model$par["theta"])) {
                    theta <- final_model$par["theta"]
                }
                
                if (!is.na(final_model$par["B0"])) {
                    B0 <- final_model$par["B0"]
                } else {
                    B0 <- NA
                }

                if (!is.na(final_model$par["lag"])) {
                    lag <- final_model$par["lag"]
                } else {
                    lag <- NA
                }

                pred <- growth_curve(final_model$par, norm_data, lag)

                print(theta)

                # Store results in the data frame
                df <- rbind(df, data.frame(
                    "raw_age" = raw_age, "fit_method" = fit_method,
                    "theta" = theta,
                    "k0" = final_model$par["k0"],
                    "B0" = B0,
                    "lag" = lag,
                    "R2" = calc_r2(norm_data, pred)
                ))
                print(paste("Iteration:", i, "raw_age:", raw_age, "fit_method:", fit_method, "theta:", theta))
            }
        }
    }
}
