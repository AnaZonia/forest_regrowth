# Gets the future predictions for different areas regrown

# select random 20% of pastureland in the Amazon
# select all pasturelands in protected areas
# select all pasturelands in the Amazon
# select all secondary forests in the Amazon

source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_cross_validate.r")
source("3_modelling/2_feature_selection.r")

# Load required libraries
library(tidyverse)
library(terra)
set.seed(1)
biome <- 1 # Amazon

model_lag <- readRDS("./0_results/amazon_model_lag.rds")
model_intercept <- readRDS("./0_results/amazon_model_intercept.rds")

# Apply Min-Max scaling using the precomputed min and max
train_stats <- readRDS("./0_results/grid_10k_amazon_secondary_train_stats.rds")

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

# pasture_selection can be random, protected, or top_20_percent

predict_future_biomass <- function(name, model, age_offset, pasture_selection = "random") {
    # Import Secondary Forest Data
    data <- import_data(paste0("grid_10k_amazon_", name), biome = biome, n_samples = "all")
    coords <- data$coords

    data <- apply_min_max_scaling(data$df, train_stats)

    if (name == "secondary") {
        data <- data %>% mutate(age = age + age_offset + round(model$par["lag"]))
    } else if (name == "pastureland") {

        # pastureland does not have an age column, so we create one
        # assuming it starts regrowing at 2020
        data$age <- age_offset
    }

    pred <- growth_curve(model$par, data, model$par["lag"])

    # get the column in data with the word "area"
    # each 1000x1000 m pixel is 1 million m2.
    # convert to hectares
    # 1 hectare = 10,000 m2
    # 1 million m2 = 100 hectares
    area <- data[[grep("area", names(data), value = TRUE)]] * 100 # convert to hectares

    # we are interested in estimated biomass gain (delta)
    if (name == "secondary") {
        pred <- pred - data$biomass
    }

    # convert biomass in Mg/ha to MgC/ha (assuming 50% C content)
    pred <- pred * 0.5

    # total estimate in Megatons of Carbon
    pred <- pred / 1000

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
        } else if (pasture_selection == "top_20_percent") {
            # select top 20% of pastureland by biomass
            # Sort predictions and data by descending prediction value
            order_indices <- order(pred, decreasing = TRUE)
            sorted_area <- area[order_indices]

            # Compute cumulative sum of area
            cum_area <- cumsum(sorted_area)

            # Find the number of pixels needed to reach 20% of total area
            total_area <- sum(area)
            n_needed <- which(cum_area >= 0.05 * total_area)[1]

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


pred_lag_2050_secondary <- predict_future_biomass("secondary", model_lag, 30)
pred_intercept_2050_secondary <- predict_future_biomass("secondary", model_intercept, 30)
pred_lag_2050_pastureland_top_20 <- predict_future_biomass("pastureland", model_lag, 30, "top_20_percent")
pred_lag_2050_pastureland_random <- predict_future_biomass("pastureland", model_lag, 30, "random")


points <- vect(pred_lag_2050_pastureland_top_20[[2]], geom = c("lon", "lat"), crs = "EPSG:4326")
writeVector(points, "0_results/pred_lag_2050_pastureland_top_20.shp", overwrite = TRUE)
points <- vect(pred_lag_2050_secondary[[2]], geom = c("lon", "lat"), crs = "EPSG:4326")
writeVector(points, "0_results/pred_lag_2050_secondary.shp", overwrite = TRUE)

# Data
df <- tibble::tibble(
    scenario = factor(
        rep(c("Secondary Forests", "Random 5% pastureland", "Top 5% pastureland"), each = 2),
        levels = c("Secondary Forests", "Random 5% pastureland", "Top 5% pastureland")
    ),
    source = rep(c("Secondary", "Pastureland"), 3),
    value = c(
        pred_lag_2050_secondary[[1]], 0,
        pred_lag_2050_secondary[[1]], pred_lag_2050_pastureland_random[[1]],
        pred_lag_2050_secondary[[1]], pred_lag_2050_pastureland_top_20[[1]]
    )
)

palette <- c(
    "Secondary" = "#4DAF4A",
    "Pastureland" = "#1F78B4"
)

p <- ggplot(df, aes(x = scenario, y = value, fill = source)) +
        geom_bar(stat = "identity", width = 0.7) +
        scale_fill_manual(values = palette, name = NULL) +
        scale_y_continuous(
            expand = expansion(mult = c(0.05, 0.05)),
            name = "Biomass stored by 2050 (Mg C)"
        ) +
        theme_minimal(base_size = 18) +
        theme(
            legend.position = "top",
            legend.text = element_text(size = 18),
            axis.text = element_text(size = 10, color = "black"),
            axis.title = element_text(size = 20, color = "black"),
            panel.grid = element_blank(),
            axis.line = element_line(color = "black", size = 1),
            axis.ticks = element_line(color = "black", size = 0.8),
            axis.ticks.length = unit(0.3, "cm"),
            plot.margin = margin(10, 10, 10, 10)
        ) +
        labs(x = NULL)

# Save to file
ggsave("0_results/figures/pastureland_biomass_stack.jpeg",
    plot = p,
    width = 8, height = 8, dpi = 300
)
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