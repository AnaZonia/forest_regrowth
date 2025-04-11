# Load necessary libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

# Get data for plotting
source("2_modelling/1_modelling.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/1_parameters.r")

biome <- 1
n_samples <- 10000
set.seed(1)
registerDoParallel(cores = 4) # Adjust cores as needed

conditions <- list('pars["k0"] < 0')

thetas <- c(0.5, 1, 2, NA)
raw_age_vals <- c(TRUE, FALSE)
fit_methods <- c("B0", "lag")

# Create an empty results list to store results from parallel execution
results_list <- foreach(theta = thetas, .combine = rbind, .packages = c("tidyverse")) %:%
    foreach(raw_age = raw_age_vals, .combine = rbind) %:%
    foreach(fit_method = fit_methods, .combine = rbind) %dopar% {
        # Initialize an empty data frame to store results
        df_local <- data.frame("raw_age" = NA, "fit_method" = NA, "theta" = NA, "k0" = NA, "B0" = NA, "lag" = NA, "R2" = NA)

        for (i in 1:5) {
            data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples) %>%
                select(-mean_pdsi)

            # Fit the model on the full data
            norm <- normalize_independently(data, raw_age)
            norm_data <- norm$train_data

            basic_pars <- c("k0", fit_method)

            if (is.na(theta)) {
                basic_pars <- c(basic_pars, "theta")
            }

            pars_init <- find_combination_pars(basic_pars, "data_pars" = c("num_fires", "sur_cover"), norm_data)
            final_model <- run_optim(norm_data, pars_init, conditions)

            theta_val <- if (!is.na(final_model$par["theta"])) final_model$par["theta"] else theta
            B0_val <- if (!is.na(final_model$par["B0"])) final_model$par["B0"] else NA
            lag_val <- if (!is.na(final_model$par["lag"])) final_model$par["lag"] else NA

            pred <- growth_curve(final_model$par, norm_data, lag_val)

            # Store results in the local data frame
            df_local <- rbind(df_local, data.frame(
                "raw_age" = raw_age, "fit_method" = fit_method,
                "theta" = theta_val, "k0" = final_model$par["k0"],
                "B0" = B0_val, "lag" = lag_val,
                "R2" = calc_r2(norm_data, pred)
            ))

            print(paste("Iteration:", i, "raw_age:", raw_age, "fit_method:", fit_method, "theta:", theta_val))
        }

        return(df_local)
    }




# Convert results list to a data frame
df <- as.data.frame(results_list)

write.csv(df, "0_results/thetas.csv", row.names = FALSE)



# get R2 with plots
# predictions for the site vs observation for the site




# regrowth trajectories for different scenarios - pick a couple examples across different plots
# add this to the map with regrowth potential

# how much of the plots are active regrowth

# show variance per age





# in the final model we combine both

# refit simultaneously with the full spread of the data

