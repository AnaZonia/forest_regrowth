####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Writing raster masks selecting for only secondary forest patches at 2019.
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAPBIOMAS raw data comes divided in 12 different regions of the Amazon, each with their own identifier (a location number).
# This script builds rasters from raw materials containing only regrowth events.
# This allows to subset the data to contain only pixels that show regrowth history
# The output masks are input into 1_processing_pred_var.R
####################################################################

library(sf)
library(terra) # handling spatial data
library(stringr) # when possible, change to tidyverse
#library(pbapply) #progress bar for apply family of functions
setwd("/home/aavila/forest_regrowth")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r") # sourcing functions

# EXTRACTING DATA FROM MULTIPLE REGIONS OF THE COUNTRY (larger scale)
path <- './mapbiomas/regrowth_raw'
files <- list.files(path)
#newname <- sub('none-','', files) ## making names easier to read, standardizing the names
#file.rename(file.path(path,files), file.path(path, newname)) ## renaming it.
locations <- str_sub(files, start= -25, end = -5)
locations <- unique(locations)
locations <- locations[5:10]

# create a raster stack with one layer per year, for this location.
# this stack will be merged to make the regrowth-only mask.

# here we go location by location:
  # import raw data
  # make a mask with only areas that show regrowth - this mask will be used to reduce the size of files we're handling
  # save these rasters in the ./regrowth_rasters directory

for (i in 1:length(locations)){
  location <- locations[1]
  print(location)
  files_tmp <- list.files(path = './mapbiomas/regrowth_raw', pattern=location, full.names=TRUE)   # obtain paths for all files for that location
  files_tmp <- sort(files_tmp)

  last_year_regrowth <- rast(files_tmp[length(files_tmp)])
  last_year_regrowth[last_year_regrowth!=303] <- NA # only leave behind values not including currently regrowing patches
  writeRaster(last_year_regrowth, file.path(paste0('./mapbiomas/regrowth_masks/', location, '_mask.tif'))) # save rasters with the year on the folder created for each location.
  
  last_year_mature <- rast(files_tmp[length(files_tmp)])
  last_year_mature[last_year_mature > 300] <- NA
  last_year_mature[last_year_mature < 200] <- NA # only leave behind values not including the regrowth moment
  writeRaster(last_year_mature, file.path(paste0('./mapbiomas/mature_masks/', location, '_mask.tif'))) # save rasters with the year on the folder created for each location.

  regrowth_list <- c()
  for (i in 1:length(files_tmp)){
    print('ith iteration of regrowth_list')
    regrowth_list <- c(regrowth_list, rast(files_tmp[i]))
    regrowth_list[[i]][regrowth_list[[i]]!=503] <- NA # only leave behind values not including the regrowth moment
  }

  writeRaster(regrowth_list[[i]], file.path(paste0('./mapbiomas/regrowth_rasters/', location, '_', c(1987+i), ".tif"))) # save rasters with the year on the folder created for each location.

  # obtain the raster file for all years within that location.
  # to make processing lighter, subset only the pixels that have shown regrowth history.
  # here, we are (1) making a mask registering all regrowth moments and (2) subsetting rasters based on that mask.
  # a regrowth moment is flagged with the value "503", therefore:

}
