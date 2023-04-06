####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Prepares fire history data
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Intakes:
#   forest_age raster
#   Raw fire rasters
# Outputs:
#   fire_history raster (time since last observed fire; total number of fires)
####################################################################

library(sf)
library(terra) # handling spatial data
library(tidyverse)
setwd("/home/aavila/forest_regrowth/dataframes")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  FIRE ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# EXTRACTING DATA FROM MULTIPLE REGIONS OF THE COUNTRY (larger scale)
path <- './mapbiomas/fire'
files <- list.files(path)
#newname <- sub('none-','', files) ## making names easier to read, standardizing the names
#file.rename(file.path(path,files), file.path(path, newname)) ## renaming it.
locations <- str_sub(files, start= -25, end = -5)
locations <- unique(locations)

#for (location in locations){}
location = '0000000000-0000095232'
regrowth_mask <- rast(paste0('./model_ready_rasters/', location, '_forest_age.tif'))

files_tmp <- list.files(path, pattern=location, full.names=TRUE)   # obtain paths for all files for that location
files_tmp <- sort(files_tmp)
tmp_rasters <- lapply(files_tmp, rast)
fire_brick <- rast(tmp_rasters)

# correcting for small inconsistencies between the sizes
fire_brick_cropped <- crop(fire_brick, ext(regrowth_mask))
regrowth_mask_cropped <- crop(regrowth_mask, ext(fire_brick_cropped))
fire_brick_masked <- mask(fire_brick_cropped, regrowth_mask_cropped)

#################################################################################
# count total number of fires
num_fires <- sum(fire_brick_masked)

# find when was last fire observed (how many years before observed regrowth)
regrowth_instances <- which.lyr(regrowth_cleaned == 503)
regrowth_last_instance <- where.max(regrowth_instances)



writeRaster(fire, paste0('./model_ready_rasters/', location, '_fire_history.tif'))
