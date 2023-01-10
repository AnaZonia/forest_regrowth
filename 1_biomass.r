####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Downloading and processing GEDI4A data.
# Ana Avila - Dec 2022
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# downloads GEDI4A product from cloud using GEDI4R package
# converts hdf5 into netcdf files for easier handling
# makes biomass dataframe from raw product with lat, lon, UTM coords, and AGBD info 
# imports and makes dataframes from santoro and potapov raw data as well.
####################################################################

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("grimbough/rhdf5")
library(rhdf5) # for handling raw GEDI data
#remotes::install_github("VangiElia/GEDI4R")
library(GEDI4R) # for extracting raw GEDI data
library(raster) # Might not need this one
library(ncdf4)
library(tidyverse)
library(rgdal)

setwd("/home/aavila/forest_regrowth")
regrowth_mask <- raster('0000000000-0000095232_mask.tif')
coords <- c(-0.5377764, -3.2823093, -48.32644, -43.99998)


# > range(age$lon)
# [1] -48.32644 -43.99998
# > range(age$lat)
# [1] -3.2823093 -0.5377764
# Dubayah et al 2022 -> GEDI L4A Footprint Level Aboveground Biomass Density (Mg/ha)
# more instructions on https://github.com/VangiElia/GEDI4R
# a package built specially for processing GEDI4A biomass data.

####################################################################
##########              BODY                 #######################
####################################################################

#shape <- readOGR(dsn = "./amazon_biome_border", layer = "amazon_biome_border")
outdir <- ("./GEDI_raw/")
#   extent(shape)[4], extent(shape)[1], extent(shape)[3], extent(shape)[2],
# 5.269581, -16.66202, -73.98318 ;  -43.39932
# > range(GEDI_mat$lat)
# [1] -8.999994 -2.000006
# > range(GEDI_mat$lon)
# [1] -63.99999 -45.00000

GEDI_download = l4_download(
  coords[1], coords[2], coords[3], coords[4],
  outdir = outdir,
  from = "2020-01-01",
  to = "2020-12-31",
  just_path = F )

# Read in
nc_data <- c(paste0("./GEDI_amazon/", list.files("./GEDI_amazon/", pattern = '.h5')))
test_convert <- function(nc_file) tryCatch(nc_open(nc_file), error = function(e) e)
nc_data2 <- lapply(nc_data, test_convert)

df_from_nc <- function(nc_file, shape){
  # Get variable names
  nc_file <- crop(nc_file, shape)
  nmv = names(nc_file$var)
  # Which ones are agbd
  nmv.values = nmv[grepl("/agbd$", nmv)]
  # Which ones are lat/lon
  nmv.lat = nmv[grepl("/lat_lowestmode$", nmv)]
  nmv.lon = nmv[grepl("/lon_lowestmode$", nmv)]

  # Get data for the first one, might want to iterate over 1:length(nmv.values), or beams or whatever
  # -9999 seems to be their way of doing NA
  df1 = data.frame(agbd = ncvar_get(nc_file, nmv.values[1]),
              lat = ncvar_get(nc_file, nmv.lat[1]),
              lon = ncvar_get(nc_file, nmv.lon[1]))

  df2 <- subset(df1, agbd > 0)
  #df3 <- subset(df2, lat < coords[1] & lat > coords[2] & lon > coords[3] & lon < coords[4])

  return(df3)
}

df_list <- lapply(nc_data2, df_from_nc, shape)

df_unified <- bind_rows(df_list)

hist(df_unified$agbd, xlim = c(0, 400), breaks <- 5000)
#looks like it starts plateauing at 10-15 ton/ha


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  SANTORO ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Santoro et al 2018 data -> GlobBiomass ESA (Mg/ha)
# 100m resolution, 2010

biomass1 = raster("./santoro/N00W060_agb.tif")
biomass2 = raster("./santoro/N00W100_agb.tif")
biomass3 = raster("./santoro/N40W060_agb.tif")

biomass = merge(biomass1, biomass2, biomass3)

## crop and mask
r2 <- crop(biomass, extent(BRA))
r3 <- mask(r2, BRA) #save this somewhere
e <- extent(-48.32644, -43.99998, -3.2823093, -0.5377764) # for testing purposes, the coordinates of region 0000000000.0000095232
r4 <- crop(biomass1, e)

bm_test <- getValues(r4)
bm_test <- data.frame(cell = 1:length(bm_test), value = bm_test)
bm_test <- na.omit(bm_test)
bm_test[,c("x","y")] <- xyFromCell(r4, bm_test$cell)

#biomass = biomass[ymin < biomass$y & biomass$y < ymax & xmin < biomass$x & biomass$x < xmax,]

biomass = cbind(biomass, LongLatToUTM(biomass$x, biomass$y))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  POTAPOV ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Potapov et al 2020 -> GLAD Forest Canopy Height (m)
# 30m resolution, 2019

biomass = raster("Forest_height_2019_SAM.tif")

## crop and mask
r2 <- crop(biomass, extent(BRA))
r3 <- mask(r2, BRA)

#writeRaster(r3, "Forest_height_2019_Brazil.tif")

e <- extent(xmin, xmax, ymin, ymax)

biomass = raster("Forest_height_2019_Brazil.tif")
biomass_cropped <- crop(biomass,e)


#biomass_df <- as.data.frame(biomass_cropped, xy=T, na.rm = TRUE)
bm_test <- values(biomass_cropped)

bm_tst_complete <- na.omit(bm_test)

bm_test <- getValues(biomass_cropped)
bm_test <- data.frame(cell = 1:length(bm_test), value = bm_test)
bm_test <- na.omit(bm_test)
bm_test[,c("x","y")] <- xyFromCell(biomass_cropped, bm_test$cell)

biomass = cbind(biomass_with_data, LongLatToUTM(biomass_with_data$x, biomass_with_data$y))