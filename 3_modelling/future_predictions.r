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

# Apply Min-Max scaling using the precomputed min and max
train_stats <- readRDS("./0_results/grid_1k_amazon_secondary_train_stats.rds")

apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Biomass predictions for future years
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Estimated biomass x years after 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# pasture_selection can be random, protected, or top_20_percent

predict_future_biomass <- function(name, model, age_offset, pasture_selection = "random") {
    # Import Secondary Forest Data
    data <- import_data(paste0("grid_1k_amazon_", name), biome = biome, n_samples = "all")
    # rename distance_deep_forest to dist
    data$df <- data$df %>% rename(dist = distance_deep_forest)
    coords <- data$coords
    data <- apply_min_max_scaling(data$df, train_stats)

    if (name == "secondary") {
        data <- data %>% mutate(age = age + age_offset)
    } else if (name == "pastureland") {
        # pastureland does not have an age column, so we create one
        # assuming it starts regrowing at 2020
        data$age <- age_offset
    }
    pred <- growth_curve(model$par, data)

    # get the column in data with the word "area"
    # each 1000x1000 m pixel is 1 million m2.
    # convert to hectares
    # 1 hectare = 10,000 m2
    # 1 million m2 = 100 hectares
    area <- data[[grep("area", names(data), value = TRUE)]] * 100 # convert to hectares

    # convert biomass in Mg/ha to C/ha (assuming 50% C content)
    pred <- pred * 0.5
    coords$pred <- pred

    if (name == "pastureland") {
        if (pasture_selection == "random") {
            set.seed(1)
            # shuffle the indices of pred
            random_indices <- sample(1:length(pred), size = 0.2 * length(pred), replace = FALSE)


            # Compute cumulative sum of area
            cum_area <- cumsum(sorted_area)

            # Find the number of pixels needed to reach 20% of total area
            total_area <- sum(area)


            pred <- pred[random_indices]
            area <- area[random_indices]
            coords <- coords[random_indices, ]
        } else if (pasture_selection == "protected") {
            # select all pasturelands in protected areas
            # get the indices where data$protedted is 1
            protected_indices <- which(data$protec == 1)
            protected_areas <- data[protected_indices, ]
            # filter coords and pred to only include protected areas
            pred <- pred[protected_indices]
            area <- area[protected_indices]
            coords <- coords[protected_indices, ]
        } else if (pasture_selection == "top_20_percent") {
            # select top 20% of pastureland by biomass
            # Sort predictions and data by descending prediction value
            order_indices <- order(pred, decreasing = TRUE)
            sorted_area <- area[order_indices]

            # Compute cumulative sum of area
            cum_area <- cumsum(sorted_area)

            # Find the number of pixels needed to reach 20% of total area
            total_area <- sum(area)
            n_needed <- which(cum_area >= 0.2 * total_area)[1]

            # Select those indices
            selected_indices <- order_indices[1:n_needed]
            pred <- pred[selected_indices]
            area <- area[selected_indices]
            coords <- coords[selected_indices, ]
        }
    }

    total_biomass <- sum(pred * area, na.rm = TRUE)

    
    return(list(total_biomass, coords))
}

pred <- growth_curve(model$par, data)
area <- data[[grep("area", names(data), value = TRUE)]] * 100 # convert to hectares



# While the sum of the selected area is less than 20% of the total area
while (sum(area) < 0.2 * sum(data$area)) {
    # Get the top 2% of predictions not already selected
    top_indices <- setdiff(top_indices[1:round(0.02 * length(pred))], selected_indices)

    # Add the corresponding areas to the total area
    area <- c(area, data$area[top_indices])

    # Keep track of the selected indices to avoid duplicates
    selected_indices <- c(selected_indices, top_indices)
}
pred <- pred[selected_indices]

# set.seed(1)
# random_indices <- sample(1:length(pred), size = 0.2 * length(pred), replace = FALSE)
# print(length(random_indices))
# pred <- pred[random_indices]
# area <- area[random_indices]
# top_20_indices <- order(pred, decreasing = TRUE)[1:round(0.2 * length(pred))]
# print(length(top_20_indices))
# pred <- pred[top_20_indices]
# area <- area[top_20_indices]
mean(area)
hist(area)



pred_lag_2050_secondary <- predict_future_biomass("secondary", model_lag, 30)
pred_intercept_2050_secondary <- predict_future_biomass("secondary", model_intercept, 30)

pred_lag_2050_pastureland_random <- predict_future_biomass("pastureland", model_lag, 30, "random")
pred_lag_2050_pastureland_protected <- predict_future_biomass("pastureland", model_lag, 30, "protected")
pred_lag_2050_pastureland_top_20 <- predict_future_biomass("pastureland", model_lag, 30, "top_20_percent")

barplot(
    c(
        # pred_lag_2050_secondary[[1]],
        # pred_intercept_2050_secondary[[1]],
        pred_lag_2050_pastureland_random[[1]],
        pred_lag_2050_pastureland_protected[[1]],
        pred_lag_2050_pastureland_top_20[[1]]
    ),
    names.arg = c(
        # "Lag 2050 Secondary",
        # "Intercept 2050 Secondary",
        "Lag 2050 Pastureland Random 20%",
        "Lag 2050 Pastureland Protected",
        "Lag 2050 Pastureland Top 20%"
    ),
    main = "Total Biomass in the Amazon (Mg C)",
    ylab = "Total Biomass (Mg C)"
)


# points <- vect(pred_lag_2050[[2]], geom = c("lon", "lat"), crs = "EPSG:4326")
# writeVector(points, "pred_lag_2050_secondary.shp", overwrite = TRUE)



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

# is more error found in slower areas?

lm(coords$percentage ~ coords$error)

# is more error found in areas of lower initial biomass?
lm(data$biomass ~ coords$error)