####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the district of Paragominas, in the Western Amazon.
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Intakes:
#   Raster mask for regrowth history
#   Raw fire rasters
# Outputs:
#   lulc dataframe (full history for all pixels showing yes/no burning detected)
#   lulc_history dataframe (time since last observed land use; total number of years under each land use)
####################################################################
library(raster)
library(data.table) #for faster reading of csv files with function fread
library(sf)
library(terra) # handling spatial data
library(tidyverse)
#library(pbapply) #progress bar for apply family of functions
setwd("/home/aavila/forest_regrowth")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")

location <- '0000000000-0000095232'

# take in mask with only locations of current regrowth
regrowth_mask <- rast(paste0('./mapbiomas/regrowth_masks/', location, '_regrowth_mask.tif'))

files <- list.files(path = './mapbiomas/lulc', pattern='\\.tif$', full.names=TRUE)   # obtain paths for all files 
tmp_rasters <- lapply(files, rast)
lulc_brick <- rast(tmp_rasters)

lulc_brick_cropped <- crop(lulc_brick, extent(regrowth_mask))
lulc_brick_masked <- mask(lulc_brick_cropped, regrowth_mask)

#################################################################################

########## LAND USE ##########
# VARIABLES OBTAINED
# number of years under each land use tnrype
# time since last observation of each land use type

# INDEX ## 3 = forest
# 15 = pasture
# 39 = soy
# 46 = coffee
# 20 = sugar cane
# 41 = other annual crop
# 48 = other perennial crop

tst <- rast('0000000000-0000095232_lulc_2.tif')
# total years under each land use type

calc_total_yrs <- function(raw_raster, val){
  raw_raster[raw_raster != val] <- 0
  total_past_years <- app(raw_raster, fun=sum)
  return(total_past_years/val)
}

pasture_total <- calc_total_yrs(tst, 15)
soy_total <- calc_total_yrs(tst, 39)
coffee_total <- calc_total_yrs(tst, 46)
sugar_total <- calc_total_yrs(tst, 20)
other_annual_total <- calc_total_yrs(tst, 41)
other_perennial_total <- calc_total_yrs(tst, 48)

# not necessary to create a function for last observed;
# just use regrowth mask and ages, and look for the land use in that position.
# we reach the column by tst[[lyr_num]][age]

pasture_instances <- which.lyr(tst == 15)
pasture_last_instance <- where.max(pasture_instances)
time_since_pasture <- max(pasture_last_instance) - pasture_last_instance

lulc <- lulc[,c(1:2, 38:ncol(lulc))]
# the time since last observed land use type, when not observed, shows as "2019" rather than NA. fixing that:
for (i in 10:12){
  lulc[,c(i)][lulc[,i] == 2019] <- NA
}

writeRaster(lulc, file.path(paste0(location, '_lulc_history.tif')))


