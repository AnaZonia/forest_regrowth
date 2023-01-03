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

library(sf)
library(raster) #  handling spatial data
library(terra) # handling spatial data
library(geodata) # to extract worldclim with getData
library(sp) # to extract worldclim with getData
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("grimbough/rhdf5")
library(rhdf5) # for handling raw GEDI data
#remotes::install_github("VangiElia/GEDI4R")
library(GEDI4R) # for extracting raw GEDI data
#library(stringr)
library(tidyverse)
library(plyr)
library(foreach) # for splitting heavy processing (masking, converting)
library(doParallel) # for splitting heavy processing (masking, converting)
## Brazil shapefile mask
library(maptools)  ## For wrld_simpl
library(​data.table​) #for faster reading of csv files with function fread
data(wrld_simpl)
BRA <- subset(wrld_simpl, NAME=="Brazil")

setwd("/home/aavila/forest_regrowth/dataframes")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")

location <- '0000000000-0000095232'
regrowth_mask <- raster(paste0('./mapbiomas/regrowth_masks/', location, '_mask.tif'))

files <- list.files(path = './mapbiomas/lulc', pattern='\\.tif$', full.names=TRUE)   # obtain paths for all files 
tmp_rasters <- lapply(files, raster)

regrowth_mask <- raster('./mapbiomas/regrowth_masks/0000000000-0000095232_mask.tif')

tmp_rasters <- pbapply::pblapply(tmp_rasters, terra::crop, extent(regrowth_mask)) # subsects all rasters to area of interest
tmp_rasters2 <-  pbapply::pblapply(tmp_rasters, terra::mask, regrowth_mask) # mask all raw files, leaving behind only the cells containing any regrowth instance.
stacked_history <- stack(tmp_rasters2)
writeRaster(stacked_tst, "0000000000-0000095232_lulc.tif")

# since, for some reason, lulc turned out with years as different bands rather than different layers (check why once Janus access is regained)
for (i in 1:35){ #total number of years
  print(i)
  lulc_raster = raster('0000000000-0000095232_lulc.tif', band = i)
  convert_history <- getValues(lulc_raster)
  convert_history <- data.frame(cell = 1:length(convert_history), value = convert_history)
  convert_history <- na.omit(convert_history)
  convert_history[,c("x","y")] <- xyFromCell(lulc_raster, convert_history$cell)
  saveRDS(convert_history, file.path(paste0('0000000000-0000095232_lulc_', c(1984+i),'.rds')))
}

files_tmp <- list.files(path = './lulc/', full.names=TRUE)   # obtain paths for all files for that location
lulc = pbapply::pblapply(files_tmp, readRDS)
lulc_df <- df_merge(lulc)

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

# total years under each land use type
lulc$pasture <- rowSums(lulc == 15)
lulc$soy <- rowSums(lulc == 39)
lulc$coffee <- rowSums(lulc == 46)
lulc$sugar <- rowSums(lulc == 20)
lulc$other_perennial <- rowSums(lulc == 48)
lulc$other_annual <- rowSums(lulc == 41)

# time since last observation of each land use type
ts_pasture <- find_last_instance(lulc, function(x) which(x == 15))
lulc$ts_pasture <- max(ts_pasture)-ts_pasture
ts_soy <- find_last_instance(lulc, function(x) which(x == 39))
lulc$ts_soy <- max(ts_soy)-ts_soy
#ts_coffee = find_last_instance(lulc, function(x) which(x == 46))
#ts_sugar = find_last_instance(lulc, function(x) which(x == 20))
ts_other_perennial <- find_last_instance(lulc, function(x) which(x == 48))
lulc$ts_other_perennial = max(ts_other_perennial)-ts_other_perennial
ts_other_annual <- find_last_instance(lulc, function(x) which(x == 41))
lulc$ts_other_annual <- max(ts_other_annual)-ts_other_annual

lulc <- cbind(lulc, LongLatToUTM(lulc$lon, lulc$lat)) # add UTM coordinates as new columns
lulc$xy <- paste0(lulc$zone, lulc$x, lulc$y)
saveRDS(fire, file.path(paste0('./mapbiomas/dataframes/', location, '_lulc_history.rds')))
