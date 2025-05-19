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
n_samples <- 10000




# Import Pastureland data

data <- import_data("unified_fc", biome = biome, n_samples = "all")
coords <- data$coords
data <- data$df

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Biomass predictions for future years
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Estimated biomass x years after 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

predict_future_biomass <- function(df, model, age_offset) {
    df <- df %>% mutate(age = age + age_offset)
    growth_curve(model$par, df, lag = model$par["lag"])
}

pred_2050 <- predict_future_biomass(norm_data, model, 30)
pred_2075 <- predict_future_biomass(norm_data, model, 55)



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

same_ages <- norm_data %>%
    mutate(age = 25)

pred_relative <- growth_curve(model$par, same_ages, lag = model$par["lag"])

coords$percentage <- coords$pred_relative / norm_data$nearest_biomass
coords$pred_relative <- pred_relative

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Percent error map for 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Find what is the estimated biomass
# when hypothetically all pixels
# are at the same age (25 years)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])

coords$percent_error <- (pred - norm_data$nearest_biomass) / norm_data$nearest_biomass





points <- vect(coords, geom = c("lon", "lat"), crs = "EPSG:4326")

# save points as a shapefile
writeVector(points, "predictions.shp", overwrite = TRUE)
