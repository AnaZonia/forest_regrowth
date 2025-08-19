# Gets the future predictions for different areas regrown

# select random 20% of pastureland in the Amazon
# select all pasturelands in protected areas
# select all pasturelands in the Amazon
# select all secondary forests in the Amazon

source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

library(tidyverse)
library(terra)
library(scales) # for label formatting

set.seed(1)
biome <- 1 # Amazon

model_lag <- readRDS("./0_results/amazon_model_lag.rds")
model_intercept <- readRDS("./0_results/amazon_model_intercept.rds")

# Apply Min-Max scaling using the precomputed min and max
train_stats <- readRDS("./0_results/grid_1k_amazon_secondary_train_stats.rds")

# keep only the variables that are in the model
train_stats <- train_stats %>%
    filter(variable %in% c(names(model_lag$par), names(model_intercept$par)))

apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}


# Biomass predictions for future years
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Estimated biomass x years after 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# pasture_selection can be random, protected, top_5_percent, or all

predict_future_biomass <- function(name, model, age_offset, pasture_selection = "random", delta = TRUE) {

    # Import Secondary Forest Data
    data <- import_data(paste0("grid_1k_amazon_", name), biome_num = biome, n_samples = "all")
    coords <- data$coords

    data <- apply_min_max_scaling(data$df, train_stats)

    data_2020 <- data

    if (name == "secondary") {
        data <- data %>% mutate(age = age + 30)
    } else if (name == "pastureland") {
        # pastureland does not have an age column, so we create one
        # assuming it starts regrowing at 2020
        data$age <- 30
        data_2020$age <- 1 # pastureland is 1 year old in 2020
    }

    pred <- growth_curve(model$par, data, model$par["lag"])

    if (delta) {
        pred_2020 <- growth_curve(model$par, data_2020, model$par["lag"])
        # if delta is TRUE, we return the difference between the two predictions
        pred <- pred - pred_2020
    }
    

    # get the column in data with the word "area"
    # each 1000x1000 m pixel is 1 million m2.
    # convert to hectares
    # 1 hectare = 10,000 m2
    # 1 million m2 = 100 hectares
    area <- data[[grep("area", names(data), value = TRUE)]] * 100 # convert to hectares

    # convert biomass in Mg/ha to MgC/ha (assuming 50% C content)
    pred <- pred * 0.5

    # total estimate in Teragram of Carbon (TgC / ha)
    pred <- pred / 1000000

    coords$pred <- pred

    if (name == "pastureland") {
        if (pasture_selection == "random") {
            set.seed(1)
            # shuffle the indices of pred
            random_indices <- sample(1:length(pred), size = 0.05 * length(pred), replace = FALSE)
            sorted_area <- area[random_indices]

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
        } else if (pasture_selection == "top_5_percent") {
            # select top 20% of pastureland by biomass
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


pred_lag_2050_secondary <- predict_future_biomass("secondary", model_lag, 30)
pred_lag_2050_pastureland_top_5 <- predict_future_biomass("pastureland", model_lag, 30, "top_5_percent")
pred_lag_2050_pastureland_random <- predict_future_biomass("pastureland", model_lag, 30, "random")
# pred_lag_2050_pastureland_all <- predict_future_biomass("pastureland", model_lag, 30, "all")

pred_lag_2050_secondary <- predict_future_biomass("secondary", model_lag, 30, delta = FALSE)
pred_lag_2050_pastureland_top_5 <- predict_future_biomass("pastureland", model_lag, 30, "top_5_percent", delta = FALSE)
pred_lag_2050_pastureland_random <- predict_future_biomass("pastureland", model_lag, 30, "random", delta = FALSE)
# pred_lag_2050_pastureland_all <- predict_future_biomass("pastureland", model_lag, 30, "all")







df <- data.frame(
    category = factor(
        c("Secondary\nForests", "Random 5% of\nPasture Cover", "Top 5% Priority\nPasture Cover"),
        levels = c("Secondary\nForests", "Random 5% of\nPasture Cover", "Top 5% Priority\nPasture Cover")
    ),
    value = c(
        pred_lag_2050_secondary[[1]],
        pred_lag_2050_pastureland_random[[1]],
        pred_lag_2050_pastureland_top_5[[1]]
    )
)


p <- ggplot(df, aes(x = category, y = value)) + # removed fill aesthetic
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
ggsave("0_results/figures/figure_4_pasture_secondary_biomass.jpeg",
    plot = p,
    width = 10, height = 12, dpi = 300
)







# make small barplot with areas

df <- data.frame(
    category = factor(
        c("Secondary\nForests", "5% of\nPasture Cover"),
        levels = c("Secondary\nForests", "5% of\nPasture Cover")
    ),
    value = c(
        pred_lag_2050_secondary[[3]],
        pred_lag_2050_pastureland_random[[3]]
    )
)

p <- ggplot(df, aes(x = category, y = value)) + # removed fill aesthetic
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
ggsave("0_results/figures/figure_4_pasture_secondary_area.jpeg",
    plot = p,
    width = 10, height = 4, dpi = 300
)



# change label top 5% priority pastureland
# keep colors black
# add legend for maps


predict_future_biomass <- function(name, model, age_offset, pasture_selection = "random", delta = TRUE) {

    # Import Secondary Forest Data
    data <- import_data(paste0("grid_1k_amazon_", name), biome_num = biome, n_samples = "all")
    coords <- data$coords

    data <- apply_min_max_scaling(data$df, train_stats)

    data_2020 <- data
    if (name == "secondary") {
        data <- data %>% mutate(age = age + 30)
    } else if (name == "pastureland") {
        data_2020$age <- 1
        data$age <- 30
    }

    pred_2050 <- growth_curve(model$par, data, model$par["lag"])

    if (delta) {
        pred_2020 <- growth_curve(model$par, data_2020, model$par["lag"])
        # if delta is TRUE, we return the difference between the two predictions
        pred <- pred - pred_2020
    }


    # get the column in data with the word "area"
    coords$pred <- pred_2050

    return(coords)
}

pred_lag_2050_pastureland_all <- predict_future_biomass("pastureland", model_lag, 30, "all")
pred_lag_2050_secondary <- predict_future_biomass("secondary", model_lag, 30)


points <- vect(pred_lag_2050_pastureland_all[[1,2]], geom = c("lon", "lat"), crs = "EPSG:4326")
writeVector(points, "0_results/pred_lag_2050_pastureland_all.shp", overwrite = TRUE)
points <- vect(pred_lag_2050_secondary[[1, 2]], geom = c("lon", "lat"), crs = "EPSG:4326")
writeVector(points, "0_results/pred_lag_2050_secondary.shp", overwrite = TRUE)



pred_lag_2050_pastureland_all <- predict_future_biomass("pastureland", model_lag, 30, "all", delta = FALSE)
pred_lag_2050_secondary <- predict_future_biomass("secondary", model_lag, 30, delta = FALSE)

points <- vect(pred_lag_2050_pastureland_all, geom = c("lon", "lat"), crs = "EPSG:4326")
writeVector(points, "0_results/total_pred_lag_2050_pastureland_all.shp", overwrite = TRUE)
points <- vect(pred_lag_2050_secondary, geom = c("lon", "lat"), crs = "EPSG:4326")
writeVector(points, "0_results/total_pred_lag_2050_secondary.shp", overwrite = TRUE)


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