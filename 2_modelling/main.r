# main.R - Main script that runs experiments

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cluster)

# Source other scripts
source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_normalize_cross_validate.r")
source("2_modelling/2_feature_selection_ga.R")
source("2_modelling/2_perm_importance.r")

# Get configuration

# Set up parallel processing
set.seed(1)
ncore = 4
registerDoParallel(cores = ncore)

# Load data
biome = 1
n_samples = 10000

data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)

# Fit the model on the full data
norm_data <- normalize_independently(data)$train_data

# Function to run a single experiment
run_experiment <- function(basic_pars_name, data_pars_name, biome) {

    # Get parameters
    basic_pars <- basic_pars_options[[basic_pars_name]]
    data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

    # Run cross-validation
    cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

    # Return summary
    return(data.frame(
        basic_pars_name = basic_pars_name,
        data_pars_name = data_pars_name,
        biome = biome,
        mean_r2 = mean(cv_results)
    ))
}

results <- data.frame()




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------------- Check Importance of parameters included
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# biome = Amazon
# basic_pars = intercept
# just vary data_pars
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

for (name in names(data_pars_options(colnames(data)))) {
    print(data_pars_options(colnames(data))[[name]])
    print("------------------------------------------------")
    result <- run_experiment("intercept", name, 1)
    print(result)
    results <- rbind(result)
}

write.csv(results, file = "results.csv", row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------------- Check Model Forms
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# biome = Amazon
# data_pars = all_mean_clim
# just vary basic_pars
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# for (name in names(basic_pars_options)) {
#     result <- run_experiment(name, "all_mean_climate", 1)
#     print(result)
#     results <- rbind(result)
# }
# n_samples <- 10000
# data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)
# # Fit the model on the full data
# norm_data <- normalize_independently(data)$train_data

# pars_init <- find_combination_pars(basic_pars = basic_pars_options[["lag"]], "data_pars" = c("num_fires", "dist"), norm_data)

# final_model <- run_optim(norm_data, pars_init, conditions)
# final_model

# pred <- growth_curve(final_model$par, data = norm_data, lag = final_model$par[["lag"]])
# r2 <- calc_r2(norm_data, pred)
# r2

# norm_data$diff <- norm_data$biomass - pred
# # plot the difference against age
# plot(norm_data$diff, norm_data$age, xlab = "Difference between observed and predicted AGB", ylab = "Age", main = "Difference vs Age")
# summary(lm(norm_data$diff ~ norm_data$age))


# tst <- unified_field_date %>%
#     # rename field_age to age
#     rename(age = field_age) %>%
#     # rename first to nearest_biomass
#     rename(nearest_biomass = first) %>%
#     # rename b1 to biomass
#     rename(biomass = field_biom) %>%
#     # keep only columns dist, num_fires, sur_cover, nearest_biomass, biomass and age
#     select(dist, num_fires, nearest_biomass, biomass, age) %>%
#     # remove rows with NA in any column
#     filter(complete.cases(.))

# tst <- normalize_independently(tst)$train_data

# pred_agb <- growth_curve(final_model$par, data = tst, lag = final_model$par[["lag"]])
# # get R2
# r2 <- calc_r2(tst, pred_agb)
# r2

# plot(pred_agb, tst$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")
# # add 1-1 line to the plot
# abline(0, 1, col = "red", lty = 2)

# # color the points by their age
# plot(pred_agb, tst$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB", col = tst$age)
# # make the color a gradient from blue to red
# library(ggplot2)
# ggplot(tst, aes(x = pred_agb, y = biomass, color = age)) +
#     geom_point() +
#     scale_color_gradient(low = "blue", high = "red") +
#     labs(x = "Predicted AGB", y = "Observed AGB", title = "Predicted vs Observed AGB") +
#     theme_minimal()

# # tell me if there is a relationship between age and distance between observed from predicted AGB
# # calculate the difference between observed and predicted AGB
# tst$diff <- tst$biomass - pred_agb
# # plot the difference against age
# plot(tst$diff, tst$age, xlab = "Difference between observed and predicted AGB", ylab = "Age", main = "Difference vs Age")
# summary(lm(tst$diff ~ tst$age))


# # check means here

# # how do I generate the value

# n_samples <- 10000
# data <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = n_samples)
# # Fit the model on the full data
# norm_data <- normalize_independently(data)$train_data

# pars_init <- find_combination_pars(basic_pars = basic_pars_options[["lag"]], "data_pars" = c("num_fires", "dist"), norm_data)

# final_model <- run_optim(norm_data, pars_init, conditions)
# final_model

# pred <- growth_curve(final_model$par, data = norm_data, lag = final_model$par[["lag"]])
# r2 <- calc_r2(norm_data, pred)
# r2

# norm_data$diff <- norm_data$biomass - pred
# # plot the difference against age
# plot(norm_data$diff, norm_data$age, xlab = "Difference between observed and predicted AGB", ylab = "Age", main = "Difference vs Age")
# summary(lm(norm_data$diff ~ norm_data$age))


# tst <- unified_field_date %>%
#     # rename field_age to age
#     rename(age = field_age) %>%
#     # rename first to nearest_biomass
#     rename(nearest_biomass = first) %>%
#     # rename b1 to biomass
#     rename(biomass = field_biom) %>%
#     # keep only columns dist, num_fires, sur_cover, nearest_biomass, biomass and age
#     select(dist, num_fires, nearest_biomass, biomass, age) %>%
#     # remove rows with NA in any column
#     filter(complete.cases(.))

# tst <- normalize_independently(tst)$train_data

# pred_agb <- growth_curve(final_model$par, data = tst, lag = final_model$par[["lag"]])
# # get R2
# r2 <- calc_r2(tst, pred_agb)
# r2

# plot(pred_agb, tst$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")
# # add 1-1 line to the plot
# abline(0, 1, col = "red", lty = 2)

# # color the points by their age
# plot(pred_agb, tst$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB", col = tst$age)
# # make the color a gradient from blue to red
# library(ggplot2)
# ggplot(tst, aes(x = pred_agb, y = biomass, color = age)) +
#     geom_point() +
#     scale_color_gradient(low = "blue", high = "red") +
#     labs(x = "Predicted AGB", y = "Observed AGB", title = "Predicted vs Observed AGB") +
#     theme_minimal()

# # tell me if there is a relationship between age and distance between observed from predicted AGB
# # calculate the difference between observed and predicted AGB
# tst$diff <- tst$biomass - pred_agb
# # plot the difference against age
# plot(tst$diff, tst$age, xlab = "Difference between observed and predicted AGB", ylab = "Age", main = "Difference vs Age")
# summary(lm(tst$diff ~ tst$age))


# # check means here

# # how do I generate the value