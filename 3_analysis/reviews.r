# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#      Run comparisons for answering to reviews
#
#                     Ana Avila - May 2026
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(foreach)
library(doParallel)
library(tidyverse)
library(xtable)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)







csv_files <- list.files(paste0("./0_data/grid_10k_secondary_GEDI_ESA_edges_removed"), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()

names(df)

# tst <- read.csv("./0_data/gedi_esa/age_gedi_esa.csv")

summary(lm(biomass ~ age, data = df))

summary(lm(ESA_biomass ~ age, data = df))

# result is the same with the data extracted for the same locations.



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
