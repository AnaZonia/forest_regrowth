# check field data

library(dplyr)
library(terra)
library(tidyverse)
library(foreach)
library(doParallel)

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

# Load data
biome <- 1
n_samples = 10000

# try it out with only 5 years old onwards to see if theta remains 1

data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000)

# Fit the model on the full data
norm_data <- normalize_independently(data)$train_data

pars_init <- find_combination_pars(basic_pars = basic_pars_options[["lag"]], "data_pars" = c("num_fires", "dist", "sur_cover"), norm_data)

final_model <- run_optim(norm_data, pars_init, conditions)

pred <- growth_curve(final_model$par, data = norm_data, lag = final_model$par["lag"])

# save plot as png
png("./0_results/extended/predicted_vs_observed_satellite.png", width = 800, height = 600)
plot(pred, norm_data$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB") # add 1-1 line to the plot
abline(0, 1, col = "red", lty = 2)
dev.off()



# ------------------------------------------------


# identify which ones are in the same site
field <- read.csv("./0_data/groa_field/field_predictors.csv")
field <- subset(field, biome == 1) # Filter for Amazon biome
field <- field %>%
    rename(biomass = field_biom,
    asymptote = nearest_mature) %>%
    mutate(age = floor(age + 0.5))

# remove rows with NA values in any colum
field <- field[complete.cases(field), ]

field <- normalize_independently(field)$train_data

# pick one random row per unique value of site_id
field_non_repeats <- field %>%
    group_by(site_id) %>%
    slice_sample(n = 1) %>%
    ungroup()
# there are 44 unique sites in the Amazon

# remove columns biome, lat, lon, site_id, plot_id
field_non_repeats <- field_non_repeats %>%
    select(-biome, -lat, -lon, -site_id, -plot_id)

pred_field <- growth_curve(final_model$par, data = norm_field)

png("./0_results/figures/extended/field_predictions_scatterplot.png", width = 800, height = 600)
plot(norm_field$age, norm_field$biomass, xlab = "Field Age", ylab = "Field Biomass", main = "Age vs Biomass")
points(norm_field$age, pred_agb, col = "red", pch = 19)
dev.off()

# get R2
r2 <- calc_r2(norm_field, pred_field)
r2

png("./0_results/figures/extended/predicted_vs_observed_field.png", width = 800, height = 600)
plot(pred_agb, norm_field$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")
# add 1-1 line to the plot
abline(0, 1, col = "red", lty = 2)
dev.off()


# ------------------------------------------------------
# Identify and handle plots with repeated measurements
# ------------------------------------------------------


# Identify plot IDs with more than one observation
plot_nums <- field %>%
    group_by(plot_id) %>%
    summarise(n = n()) %>%
    filter(n > 1) %>%
    arrange(desc(n))

# take only rows with more than 5 observations for the same plot_id
field_repeats <- field_repeats %>%
    group_by(plot_id) %>%
    filter(n() > 5) %>%
    ungroup()

pred_repeats
table(field_repeats$site_id)
pred_repeats <- growth_curve(final_model$par, data = field_repeats)
plot(field_repeats$age, field_repeats$biomass, xlab = "Age", ylab = "Biomass", main = "Predictions for site with repeated measurements")
points(field_repeats$age, pred_repeats, col = "red", pch = 19)

