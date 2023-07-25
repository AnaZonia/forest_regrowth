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

library(terra) # handling spatial data
library(tidyverse)
setwd("/home/aavila/forest_regrowth")

extent <- c(-47.7, -46.6, -3.4, -2.6)
fire <- rast('./model_ready_rasters/0000000000-0000095232_fire_history_santoro.tif')
# take in mask with only locations of current regrowth
#regrowth_mask <- rast(paste0('./model_ready_rasters/', location, '_forest_age.tif'))
regrowth_mask <- rast('./secondary_forest_age_v2_2018/secondary_forest_age_v2_2018-0000000000-0000065536.tif')
regrowth_mask[regrowth_mask == 0] <- NA # not secondary; regrowth_mask == 1 means conversion in 2017
regrowth_mask[regrowth_mask == 33] <- NA #conversion in 1985
regrowth_mask <- crop(regrowth_mask, fire)

files <- list.files(path = './mapbiomas/lulc', pattern='\\.tif$', full.names=TRUE)   # obtain paths for all files 
tmp_rasters <- lapply(files, rast)
tmp_rasters <- lapply(tmp_rasters, crop, regrowth_mask)
tmp_rasters <- lapply(tmp_rasters, mask, regrowth_mask)
tmp_rasters <- tmp_rasters[1:33] #remove 2019 and 2018 - 2017 conversion is a forest of age 1

zerotoNA <- function(x){
  x[x==0]<-NA
  return(x)}
tmp_rasters <- lapply(tmp_rasters, zerotoNA)

#---------------------
# if age is 1 in 2018, it was not forest in 2017
# which means I want category in 2017
# if age is 33 in 2018, it was not forest in 1985
# which means I want category in 1985

# last observed land use type
filename_indx <- paste0(tempfile(), "_.tif")
layer_indices_tiles <- makeTiles(regrowth_mask, c(2000,2000), filename_indx, na.rm = TRUE, overwrite=TRUE)

row_indices <- rep(1:2000, each = 2000) 
col_indices <- rep(1:2000, times = 2000) 

get_last_layer_val <- function(lu_rast, regrowth, val){
  layer_indices_mat <- as.matrix(layer_indices_rst, wide=TRUE)
  layer_indices <- as.vector(layer_indices_rst)

  layer_index <- layer_indices
  layer_index[layer_index != val] <- NA
  tst <- tmp_rasters[[1]][cbind(row_indices, col_indices, layer_index)]

  filename <- paste0(tempfile(), "_.tif")
  ff <- makeTiles(tmp_rasters2[[10]], c(2000,2000), filename, na.rm = TRUE, overwrite=TRUE)

  layer_index <- layer_indices
  layer_index[layer_index != val] <- NA
  tst <- tmp_rasters[[1]][cbind(row_indices, col_indices, layer_index)]

  return(lu_rast[cbind(row_indices, col_indices, layer_index)])
}


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
calc_total_yrs <- function(masked_brick, val){
  masked_brick[masked_brick != val] <- 0
  total_past_years <- sum(masked_brick)
  return(total_past_years/val)
}

pasture_total <- calc_total_yrs(lulc_brick_masked, 15)
soy_total <- calc_total_yrs(lulc_brick_masked, 39)
coffee_total <- calc_total_yrs(lulc_brick_masked, 46)
sugar_total <- calc_total_yrs(lulc_brick_masked, 20)
other_annual_total <- calc_total_yrs(lulc_brick_masked, 41)
other_perennial_total <- calc_total_yrs(lulc_brick_masked, 48)

# not necessary to create a function for last observed;
# just use regrowth mask and ages, and look for the land use in that position.
# we reach the column by tst[[lyr_num]][age]

calc_time_since_lu <- function(masked_brick, val){
  masked_brick_flipped <- masked_brick [[ c(rev(order(names(masked_brick)))) ]]
  lu_instances <- which.lyr(masked_brick == val) # position of layer with last observation of land use type designated by value val
  return(nlyr(masked_brick_flipped) - lu_instances)
}

ts_pasture <- calc_time_since_lu(lulc_brick_masked, 15)
ts_soy <- calc_time_since_lu(lulc_brick_masked, 39)
ts_coffee <- calc_time_since_lu(lulc_brick_masked, 46)
ts_sugar <- calc_time_since_lu(lulc_brick_masked, 20)
ts_other_annual <- calc_time_since_lu(lulc_brick_masked, 41)
ts_other_perennial <- calc_time_since_lu(lulc_brick_masked, 48)

# note that there are land use types in last_LU that are not the 4 most common (soy, pasture, other perennial, other annual)
# decided to keep them for now but can be removed later with this:
# last_LU[last_LU != 15 & last_LU != 39 & last_LU != 41 & last_LU != 48] <- NA

lulc <- c(last_LU, ts_pasture, ts_soy, ts_other_annual, ts_other_perennial,
pasture_total, soy_total, other_annual_total, other_perennial_total)
names(lulc) <- c('last_LU', 'ts_pasture', 'ts_soy', 'ts_other_annual', 'ts_other_perennial',
'pasture_total', 'soy_total', 'other_annual_total', 'other_perennial_total')

writeRaster(lulc, paste0('./model_ready_rasters/', location, '_lulc_history_santoro.tif'))
