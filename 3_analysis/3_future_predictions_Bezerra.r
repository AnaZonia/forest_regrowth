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
# --------- Future predictions - Bezerra map ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# make future projections with the general model from error propagation
# get the projection for each pixel selected from the 10km grid cell
# get the area expected for that grid cell, then from that get the total

# get average carbon per area for those cells

# multiply that by the total area regrown to get the estimated

# ----------------------------------------------------

library(terra)

# veg <- forest vegetation
# gveg <- grassland vegetation
# mosc <- mosaic vegetation
# fores <- forestry

for (scenario in c("SSP1_RCP19", "SSP2_RCP45", "SSP3_RCP70")) {
    r <- rast(paste0("./0_data/LUCCMEBR_", scenario, "_land_cover_type_100km2_2015_2050.nc"))

    forest_2050 <- r$veg_8
    forest_2015 <- r$veg_1

    writeRaster(forest_2050, paste0("./0_data/forest_2050_", scenario, ".tif"), overwrite = TRUE)
    writeRaster(forest_2015, paste0("./0_data/forest_2015_", scenario, ".tif"), overwrite = TRUE)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------- Calculate future carbon sequestration ---------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        print(var)
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}

error_prop_results <- readRDS("./0_results/0_error_prop.rds")

pars <- error_prop_results[[2]]
train_stats <- error_prop_results[[3]]
train_stats <- subset(train_stats, variable != "sd")


future <- read.csv("./0_data/future_scenarios_area.csv")
# remove rows where cec is NA (deal with this later)
future <- subset(future, !is.na(cec))
norm_future <- apply_min_max_scaling(future[, names(future) %in% train_stats$variable], train_stats)
future <- future %>% select(growth_SSP1_RCP19, growth_SSP2_RCP45, growth_SSP3_RCP70, nearest_mature, topography, floodable_forests, protec, indig)
future <- cbind(future, norm_future)

future$age <- 35
future <- future %>%
    rename(asymptote = nearest_mature)

future <- dummy_cols(future,
    select_columns = "topography",
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)

# the area is in km2. To convert to hectares, multiply by 100

future[, c(1:3)] <- future[, c(1:3)]*100

ssp1 <- subset(future, !is.na(growth_SSP1_RCP19))
ssp2 <- subset(future, !is.na(growth_SSP2_RCP45))
ssp3 <- subset(future, !is.na(growth_SSP3_RCP70))

ssp1$pred <- growth_curve(pars, ssp1, pars["lag"])
ssp2$pred <- growth_curve(pars, ssp2, pars["lag"])
ssp3$pred <- growth_curve(pars, ssp3, pars["lag"])

ssp1 <- subset(ssp1, !is.na(pred))
ssp2 <- subset(ssp2, !is.na(pred))
ssp3 <- subset(ssp3, !is.na(pred))

ssp1_total <- sum(ssp1$pred)
ssp2_total <- sum(ssp2$pred)
ssp3_total <- sum(ssp3$pred)

ssp1_total_area <- sum(ssp1$growth_SSP1_RCP19)
ssp2_total_area <- sum(ssp2$growth_SSP2_RCP45)
ssp3_total_area <- sum(ssp3$growth_SSP3_RCP70)



pars <- pars[names(pars)[!names(pars) %in% c("floodable_forests")]]

data_1k <- import_data("grid_1k_amazon_pastureland", biome = 1, n_samples = "all")

coords <- data_1k$coords
data_1k <- data_1k$df

data_1k <- apply_min_max_scaling(data_1k, train_stats)

data_2015 <- data_1k
# pastureland does not have an age column, so we create one
# assuming it starts regrowing at 2015
data_1k$age <- 35
data_2015$age <- 1 # pastureland is 1 year old in 2020
# pastureland data here is artificially set to 30 - it isn't lagged! adding lag here will make it artificially smaller
pred <- growth_curve(pars, data_1k)
pred_2020 <- growth_curve(pars, data_2015)


area <- data_1k[[grep("area", names(data_1k), value = TRUE)]] / 10000 # convert to hectares
# get total biomass (assuming 25.8% of total biomass is belowground)
pred <- pred + (pred * 0.258)
# convert biomass in Mg/ha to MgC/ha (assuming 50% C content)
pred <- pred * 0.5
# total estimate in Teragram of Carbon (TgC / ha)
pred <- pred / 1000000





# select top ssp1_total_area of pastureland by biomass
# Sort predictions and data by descending prediction value
order_indices <- order(pred, decreasing = TRUE)
sorted_area <- area[order_indices]

# Compute cumulative sum of area
cum_area <- cumsum(sorted_area)

# Find the number of pixels needed to reach 5% of total area
total_area <- sum(area)
n_needed <- which(cum_area >= ssp1_total_area)[1]

# Select those indices
selected_indices <- order_indices[1:n_needed]
pred <- pred[selected_indices]
area <- area[selected_indices]
coords <- coords[selected_indices, ]


