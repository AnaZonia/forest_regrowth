
# Load necessary libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

# get data for plotting
source("2_R_scripts/1_modelling.r")
source("2_R_scripts/1_data_processing.r")
source("2_R_scripts/1_parameters.r")

biome <- 1
n_samples <- 10000
set.seed(1)
registerDoParallel(cores = 4)


source("2_R_scripts/modelling_no_theta.r")

conditions <- list('pars["k0"] < 0')

print("Unnormalized age")
print("lag")
print("theta = 2")
#

for (i in 1:5) {
    data <- import_data("./0_data/unified_fc_old_biomass.csv", biome = biome, n_samples = n_samples) %>%
        rename(biomass = b1) %>%
        # remove mean_pdsi column
        select(-mean_pdsi)

    # Fit the model on the full data
    norm <- normalize_independently(data)
    norm_data <- norm$train_data
    norm_stats <- norm$train_stats
    # print(norm_stats)

    head(norm_data[, c("age", "num_fires", "sur_cover", "nearest_biomass", "biomass")])

    # use optim GA

    pars_init <- find_combination_pars(basic_pars = c("lag", "k0", "theta"), "data_pars" = c("num_fires", "sur_cover"), norm_data)

    print(pars_init)

    final_model <- run_optim(norm_data, pars_init, conditions)

    print(final_model$par)
    pred <- growth_curve(final_model$par, norm_data, final_model$par[["lag"]])
    print(calc_r2(norm_data, pred))

}


# check the effect of different theta values
# check whether age is normalized or not
# check whether that makes a difference in lag or B0 model


# for each of these, I need final R2, the initial parameters, and final parameters



