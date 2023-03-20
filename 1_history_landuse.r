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

library(data.table) #for faster reading of csv files with function fread
library(sf)
library(terra) # handling spatial data
library(tidyverse)
library(pbapply) #progress bar for apply family of functions
setwd("/home/aavila/forest_regrowth")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")

location <- '0000000000-0000095232'
regrowth_mask <- rast(paste0('./mapbiomas/regrowth_masks/', location, '_mask.tif'))

files <- list.files(path = './mapbiomas/lulc', pattern='\\.tif$', full.names=TRUE)   # obtain paths for all files 
tmp_rasters <- lapply(files, raster)

regrowth_mask <- raster('./mapbiomas/regrowth_masks/0000000000-0000095232_mask.tif')

tmp_rasters <- pbapply::pblapply(tmp_rasters, terra::crop, extent(regrowth_mask)) # subsects all rasters to area of interest
tmp_rasters2 <-  pbapply::pblapply(tmp_rasters, terra::mask, regrowth_mask) # mask all raw files, leaving behind only the cells containing any regrowth instance.
stacked_history <- brick(tmp_rasters2)
writeRaster(stacked_tst, "0000000000-0000095232_lulc.tif")

lulc_raster = raster('0000000000-0000095232_lulc.tif', band = i)
convert_history <- getValues(lulc_raster)
convert_history <- data.frame(cell = 1:length(convert_history), value = convert_history)
convert_history <- na.omit(convert_history)
convert_history[,c("x","y")] <- xyFromCell(lulc_raster, convert_history$cell)
saveRDS(convert_history, file.path(paste0('0000000000-0000095232_lulc_', c(1984+i),'.rds')))



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
lulc$pasture <- rowSums(lulc[,c(3:37)] == 15)
lulc$soy <- rowSums(lulc[,c(3:37)] == 39)
lulc$coffee <- rowSums(lulc[,c(3:37)] == 46)
lulc$sugar <- rowSums(lulc[,c(3:37)] == 20)
lulc$other_perennial <- rowSums(lulc[,c(3:37)] == 48)
lulc$other_annual <- rowSums(lulc[,c(3:37)] == 41)

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

# the column names all came out as 'unlist.last.', so we change them here
colnames(lulc[,c(9:12)]) <- c('ts_pasture', 'ts_soy', 'ts_other_perennial', 'ts_other_annual')
lulc <- lulc[,c(1:2, 38:ncol(lulc))]
# the time since last observed land use type, when not observed, shows as "2019" rather than NA. fixing that:
for (i in 10:12){
  lulc[,c(i)][lulc[,i] == 2019] <- NA
}

lulc <- cbind(lulc, LongLatToUTM(lulc$lon, lulc$lat)) # add UTM coordinates as new columns
lulc$xy <- paste0(lulc$zone, lulc$x, lulc$y)
saveRDS(lulc, file.path(paste0(location, '_lulc_history.rds')))


