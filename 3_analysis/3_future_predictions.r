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

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------- Train model with 10k dataset ------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data_10k <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000)

data_10k_scaled <- normalize_independently(data_10k)
train_stats <- data_10k_scaled$train_stats
data_10k_scaled <- data_10k_scaled$train_data

init_params <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = setdiff(data_pars_options(colnames(data_10k))[["all_mean_climate"]], "floodable_forests"),
    data = data_10k_scaled
)

model <- run_optim(data_10k_scaled, init_params, conditions)

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

predict_future_biomass <- function(name, pasture_selection = "random", age_offset = 30, delta = TRUE) {

    data_1k <- import_data(paste0("grid_1k_amazon_", name), biome_num = 1, n_samples = "all")
    coords <- data_1k$coords

    data_1k <- apply_min_max_scaling(data_1k$df, train_stats)

    data_2020 <- data_1k

    # we want the values 30 years in the future
    if (name == "secondary") {
        data_1k <- data_1k %>% mutate(age = age + age_offset)
    } else if (name == "pastureland") {
        # pastureland does not have an age column, so we create one
        # assuming it starts regrowing at 2020
        data_1k$age <- age_offset
        data_2020$age <- 1 # pastureland is 1 year old in 2020
    }

    pred <- growth_curve(model$par, data_1k, model$par["lag"])

    if (delta) {
        pred_2020 <- growth_curve(model$par, data_2020, model$par["lag"])
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

    area <- data_1k[[grep("area", names(data_1k), value = TRUE)]] * 100 # convert to hectares

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
        } else if (pasture_selection == "protected") {
            # select all pasturelands in protected areas
            # get the indices where data$protedted is 1
            protected_indices <- which(data$protec == 1)
            protected_areas <- data[protected_indices, ]
            # filter coords and pred to only include protected areas
            pred <- pred[protected_indices]
            area <- area[protected_indices]
            coords <- coords[protected_indices, ]
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


pred_2050_secondary <- predict_future_biomass("secondary")
pred_2050_pastureland_top_5 <- predict_future_biomass("pastureland", "top_5_percent")
pred_2050_pastureland_random <- predict_future_biomass("pastureland", "random")
pred_2050_pastureland_all <- predict_future_biomass("pastureland", "all")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------- Figure 4 c ----------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---- Barplot with areas of sec. forests and pastures ---- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

fig_4_c <- data.frame(
    category = factor(
        c("Secondary\nForests", "5% of\nPasture Cover"),
        levels = c("Secondary\nForests", "5% of\nPasture Cover")
    ),
    value = c(
        pred_2050_secondary[[3]],
        pred_2050_pastureland_random[[3]]
    )
)

fig_4_c <- ggplot(fig_4_c, aes(x = category, y = value)) + # removed fill aesthetic
    geom_bar(stat = "identity", width = 0.5, fill = "black") + # set fill manually
    scale_y_continuous(
        labels = scales::label_comma(),
        name = "Area covered (million hectares)"
    ) +
    labs(x = NULL) +
    theme_minimal(base_size = 8) +
    theme(
        axis.title.x = element_text(size = 25, color = "black", margin = margin(t = 15)),
        axis.text.x = element_text(size = 30, color = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 25, color = "black", margin = margin(r = 15)),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5)
    ) +
    coord_flip()


# Save to file
ggsave("0_results/figures/figure_4_c.jpeg",
    plot = fig_4_c,
    width = 10, height = 4, dpi = 300
)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------- Figure 4 d ----------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------- Barplot with carbon sequestered by 2050 -------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

fig_4_d <- data.frame(
    category = factor(
        c("Secondary\nForests",
        "Random 5% of\nPasture Cover",
        "Top 5% Priority\nPasture Cover"),
        levels = c("Secondary\nForests", 
            "Random 5% of\nPasture Cover", 
            "Top 5% Priority\nPasture Cover")
    ),
    value = c(
        pred_2050_secondary[[1]],
        pred_2050_pastureland_random[[1]],
        pred_2050_pastureland_top_5[[1]]
    )
)


fig_4_d <- ggplot(fig_4_d, aes(x = category, y = value)) + # removed fill aesthetic
    geom_bar(stat = "identity", width = 0.7, fill = "black") + # set fill manually
    scale_y_continuous(
        labels = scales::label_comma(),
        name = "Carbon stored by 2050 (Tg CO₂e)"
    ) +
    labs(x = NULL) +
    theme_minimal(base_size = 16) +
    theme(
        axis.title.y = element_text(size = 30, color = "black", margin = margin(r = 15)),
        axis.text.x = element_text(size = 25, color = "black", margin = margin(t = 15)),
        axis.text.y = element_text(size = 25, color = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.8)
    )

# Save to file
ggsave("0_results/figures/figure_4_d.jpeg",
    plot = fig_4_d,
    width = 10, height = 12, dpi = 300
)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------- Export maps ---------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

writeVector(
    vect(pred_2050_pastureland_all[[2]], geom = c("lon", "lat"), crs = "EPSG:4326"),
    "0_results/figures/QGIS/predictions/pred_2050_pastureland_all.shp",
    overwrite = TRUE
)

writeVector(
    vect(pred_2050_secondary[[2]], geom = c("lon", "lat"), crs = "EPSG:4326"), "0_results/figures/QGIS/predictions/pred_2050_secondary.shp",
    overwrite = TRUE
)









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