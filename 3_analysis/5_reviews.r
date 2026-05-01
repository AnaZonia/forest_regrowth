library(tidyverse)

csv_files <- list.files(paste0("./0_data/gedi_esa"), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()


tst <- read.csv("./0_data/gedi_esa/age_gedi_esa.csv")

summary(lm(biomass ~ age, data = tst))

summary(lm(ESA_biomass ~ age, data = tst))


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
