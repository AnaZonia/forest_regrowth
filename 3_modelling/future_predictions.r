# Gets the future predictions for different areas regrown

# select random 20% of pastureland in the Amazon
# select all pasturelands in protected areas
# select all pasturelands in the Amazon
# select all secondary forests in the Amazon

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(terra)

source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_normalize_cross_validate.r")
source("3_modelling/2_feature_selection_ga.R")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

biome <- 1

model_lag <- readRDS("./0_results/amazon_model_lag.rds")
model_intercept <- readRDS("./0_results/amazon_model_intercept.rds")

# Import Secondary Forest Data
data <- import_data("grid_1k_amazon_secondary", biome = biome, n_samples = "all")
coords <- data$coords
data <- data$df

# Apply Min-Max scaling using the precomputed min and max
train_stats <- readRDS("./0_results/grid_1k_amazon_secondary_train_stats.rds")

for (i in seq_along(train_stats$variable)) {
    var <- train_stats$variable[i]
    print(var)
    data[[var]] <- (data[[var]] - train_stats$min[i]) /
        (train_stats$max[i] - train_stats$min[i])
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Biomass predictions for future years
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Estimated biomass x years after 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

predict_future_biomass <- function(data, model, age_offset) {
    data <- data %>% mutate(age = age + age_offset)
    growth_curve(model$par, data)
}

pred_lag_2050 <- predict_future_biomass(data, model_lag, 30)
pred_lag_2075 <- predict_future_biomass(data, model_lag, 55)

pred_intercept_2050 <- predict_future_biomass(data, model_intercept, 30)
pred_intercept_2075 <- predict_future_biomass(data, model_intercept, 55)


coords$pred_lag_2050 <- pred_lag_2050
points <- vect(coords, geom = c("lon", "lat"), crs = "EPSG:4326")
writeVector(points, "pred_lag_2050.shp", overwrite = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ~~~~~~~~~~~~~~~~~~~~ Whole Amazon ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check Relative Growth Rate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Find what is the estimated biomass
# when hypothetically all pixels
# are at the same age (25 years)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# remove intercept to get relative growth rate

relative_growth_rate <- function(data, model) {
    same_ages <- data %>%
        mutate(age = 40)
    growth_curve(model$par, same_ages)
}

relative_lag <- relative_growth_rate(data, model_lag)
relative_intercept <- relative_growth_rate(data, model_intercept)


coords$percentage <- (relative_lag * 100) / data$nearest_mature
coords$pred_40_years <- relative_lag
points <- vect(coords, geom = c("lon", "lat"), crs = "EPSG:4326")
# save points as a shapefile
writeVector(points, "relative_growth.shp", overwrite = TRUE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Error map for 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Find what is the estimated biomass
# when hypothetically all pixels
# are at the same age (25 years)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# get error in magnitude
# add 1 to the data$biomass with 0 values
# data$biomass[data$biomass == 0] <- 1
data$biomass[data$biomass == 0] <- 1

pred <- growth_curve(model_lag$par, data, model_lag$par["lag"])

coords$error <- ((pred - data$biomass) / data$biomass) * 100

points <- vect(coords, geom = c("lon", "lat"), crs = "EPSG:4326")
# save points as a shapefile
writeVector(points, "percent_error.shp", overwrite = TRUE)
