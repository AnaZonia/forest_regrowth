# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#    Predictions for Regrowth by 2050 (pasture and secondary)
#
#                 Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

library(tidyverse)
library(terra)
library(scales) # for label formatting
library(foreach)
library(doParallel)

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Estimate biomass by 2050 --------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'
#' @param name Character. The land cover type: "secondary" or "pastureland".
#' @param age_offset Integer. Number of years in the future to predict (e.g., 30 for 2020 → 2050).
#' @param pasture_selection Character. How to select pasture pixels:
#'        "random": Default. Randomly selects 5% of pastureland pixels.
#'        "protected": Pasturelands in protected areas.
#'        "top_5_percent": 5% of pastureland area with highest regrowth potential.
#'        "all": All pixels.
#' @param delta Logical. If TRUE, returns only the increment in biomass over the forecast period; if FALSE, predicts total value in target year.
#'
#' @return List with components:
#'   - total_biomass: Predicted total biomass (TgC)
#'   - coords: Data frame of pixel coordinates and predictions
#'   - total_area: Total area included in the sum (million hectares)


predict_future_biomass <- function(name, model, train_stats, pasture_selection = "random", age_offset = 30, delta = TRUE) {

    data_1k <- import_data(paste0("grid_1k_amazon_", name), biome = 1, n_samples = "all")

    coords <- data_1k$coords
    data_1k <- data_1k$df

    data_1k <- apply_min_max_scaling(data_1k, train_stats)

    data_2020 <- data_1k

    # we want the values 30 years in the future
    if (name == "secondary") {
        data_1k <- data_1k %>% mutate(age = age + age_offset)
        # include the lag term here since we aren't trying to find the biomass for those pixels if they were starting from zero.
        # they start from the lag-corrected age (25 + age), so that's the same data that would be used for R2 estimation
        # if they started from zero, the total biomass would be lesser, but the delta from zero to max would have been greater
        pred <- growth_curve(model$par, data_1k, model$par["lag"])
        pred_2020 <- growth_curve(model$par, data_2020, model$par["lag"])
    } else if (name == "pastureland") {
        # pastureland does not have an age column, so we create one
        # assuming it starts regrowing at 2020
        data_1k$age <- age_offset
        data_2020$age <- 1 # pastureland is 1 year old in 2020
        # pastureland data here is artificially set to 30 - it isn't lagged! adding lag here will make it artificially smaller
        pred <- growth_curve(model$par, data_1k)
        pred_2020 <- growth_curve(model$par, data_2020)
    }

    if (delta) {
        # if delta is TRUE, we return the difference between the two predictions
        # otherwise, we return the prediction for 2050
        pred <- pred - pred_2020
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # ----- Get the "area" column and convert to hectares ----- #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1000x1000 m pixel = 1 million m2 = 100 hectares
    # (1 hectare = 10,000 m2)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    area <- data_1k[[grep("area", names(data_1k), value = TRUE)]] / 10000 # convert to hectares

    # get total biomass (assuming 25.8% of total biomass is belowground)
    pred = pred + (pred * 0.258)

    # convert biomass in Mg/ha to MgC/ha (assuming 50% C content)
    pred <- pred * 0.5

    # total estimate in Teragram of Carbon (TgC / ha)
    pred <- pred / 1000000

    coords$pred <- pred

    if (name == "pastureland") {

        if (pasture_selection == "random") {
            # shuffle the indices of pred
            random_indices <- sample(1:length(pred), size = 0.05 * length(pred), replace = FALSE)

            sorted_area <- area[random_indices]
            # Compute cumulative sum of area
            cum_area <- cumsum(sorted_area)
            # Find the number of pixels needed to reach 5% of total area
            total_area <- sum(area)

            pred <- pred[random_indices]
            area <- area[random_indices]
            coords <- coords[random_indices, ]

        } else if (pasture_selection == "top_5_percent") {
            # select top 5% of pastureland by biomass
            # Sort predictions and data by descending prediction value
            order_indices <- order(pred, decreasing = TRUE)
            sorted_area <- area[order_indices]

            # Compute cumulative sum of area
            cum_area <- cumsum(sorted_area)

            # Find the number of pixels needed to reach 5% of total area
            total_area <- sum(area)
            n_needed <- which(cum_area >= 0.05 * total_area)[1]

            # Select those indices
            selected_indices <- order_indices[1:n_needed]
            pred <- pred[selected_indices]
            area <- area[selected_indices]
            coords <- coords[selected_indices, ]
        } else if (pasture_selection == "all") {
            # use all pastureland
            # nothing to do here, pred, area, and coords are already set
        }
    }

    total_biomass <- sum(pred * area, na.rm = TRUE)
    total_area <- sum(area, na.rm = TRUE) / 1000000 # convert to million hectares

    return(list(total_biomass, coords, total_area))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------- Train model with 10k dataset ------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 30000)
indices <- sample(c(1:5), nrow(data), replace = TRUE)

pred_2050_secondary_list <- numeric(5)
pred_2050_secondary_df <- data.frame()

scenarios <- c("top_5_percent", "random", "all")
carbon_lists <- list()
area_lists <- list()

for (scenario in scenarios) {
    carbon_lists[[scenario]] <- numeric(5)
    area_lists[[scenario]] <- numeric(5)
}

pred_2050_pastureland_all_df <- data.frame()


for (index in 1:5) {
    train_data <- data[indices == index, ]

    norm_data <- normalize_independently(train_data)
    train_stats <- norm_data$train_stats
    norm_data <- norm_data$train_data

    pars_init <- find_combination_pars(
        basic_pars = basic_pars_options[["lag"]],
        data_pars = setdiff(data_pars_options(colnames(data))[["all"]], "floodable_forests"),
        data = norm_data
    )

    model <- run_optim(norm_data, pars_init[[1]], conditions)

    pred_2050_secondary <- predict_future_biomass("secondary", model, train_stats)
    pred_2050_secondary_list[index] <- pred_2050_secondary[[1]]

    for (scenario in scenarios) {
        prediction <- predict_future_biomass("pastureland", model, train_stats, scenario)
        carbon_lists[[scenario]][index] <- prediction[[1]]
        area_lists[[scenario]][index] <- prediction[[3]]

        if (scenario == "all") {
            if (index == 1) {
                pred_2050_pastureland_all_df <- prediction[[2]]
                pred_2050_secondary_df <- pred_2050_secondary[[2]]
            } else {
                pred_2050_pastureland_all_df <- cbind(pred_2050_pastureland_all_df, prediction[[2]]$pred)
                pred_2050_secondary_df <- cbind(pred_2050_secondary_df, pred_2050_secondary[[2]]$pred)
            }
        }
    }
    
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ---------------------- Export maps ----------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Helper function to:
# 1) Compute mean of columns 3:7
# 2) Convert Tg/ha to Mg/ha (x 1e6)
# 3) Keep only lon, lat, pred

prepare_prediction_df <- function(df) {
    df$pred <- rowMeans(df[, 3:7], na.rm = TRUE)
    df$pred <- df$pred * 1e6 # Tg/ha -> Mg/ha
    df <- df[, c("lon", "lat", "pred")]
    return(df)
}

# Prepare pastureland predictions
pred_2050_pastureland_all_df <- prepare_prediction_df(pred_2050_pastureland_all_df)

# Prepare secondary vegetation predictions
pred_2050_secondary_df <- prepare_prediction_df(pred_2050_secondary_df)



range(pred_2050_pastureland_all_df$pred)




# Export to shapefiles
writeVector(
    vect(pred_2050_pastureland_all_df, geom = c("lon", "lat"), crs = "EPSG:4326"),
    "0_results/figures/QGIS/predictions/pred_2050_pastureland_all.shp",
    overwrite = TRUE
)
