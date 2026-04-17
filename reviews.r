

library(tidyverse)

csv_files <- list.files(paste0("./0_data/grid_10k_amazon_secondary_allpixels"), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()

head(df)

unique(df$ecoreg)

table(df$edge)

# df <- df %>%
#     filter(ecoreg == 518)

# nrow(df)

# edge_pixels <- df %>%
#     filter(edge == 1)

# non_edge_pixels <- df %>%
#     filter(edge == 0)

# mean(edge_pixels$biomass, na.rm = TRUE)
# mean(non_edge_pixels$biomass, na.rm = TRUE)


csv_files <- list.files(paste0("./0_data/esacci_sd"), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()

head(df)

plot(df$sd, df$biomass)




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


