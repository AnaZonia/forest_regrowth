####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Extracting regrowth information for the Brazilian Amazon
# Ana Avila - Jan 2023
# Intakes:
#   Raster mask with currently regrowing forests
#   Raw regrowth rasters
# Outputs:
#   forest_age raster (time of last observed regrowth)
####################################################################

library(sf)
library(terra) # handling spatial data
library(tidyverse)
setwd("/home/aavila/forest_regrowth")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r") # sourcing functions
location = '0000000000-0000095232'

#for (location in locations){}
files_tmp <- list.files(path = './mapbiomas/regrowth', pattern=location, full.names=TRUE)   # obtain paths for all files for that location
files_tmp <- sort(files_tmp)
regrowth_mask <- rast(paste0('./mapbiomas/regrowth_masks/', location, '_regrowth_mask.tif'))
# 3.002014 % are classified as currently regrowing forests with the code 303 (which seems to be the most common)

tmp_rasters <- lapply(files_tmp, rast)
reg_brick <- rast(tmp_rasters) # makes it into a brick

reg_brick_masked <- mask(reg_brick_cropped, regrowth_mask)

# now, we use the regrowth mask to make a regrowth stack with all history per location,
#containing data only on pixels that show any regrowth event in their history.

###################################
########## DATA CLEANING ##########
###################################

# INDEX #
# 100 anthropic
# 200 mature
# 300 secondary
# 400 deforestation
# 500 regrowth
# 600 secondary suppression

reg_brick_masked[reg_brick_masked < 400 & reg_brick_masked >= 300] <- 0 #secondary
reg_brick_masked[100 <= reg_brick_masked & reg_brick_masked < 200] <- 1 #anthropic
reg_brick_masked[reg_brick_masked == 515] <- 1 #anthropic - fixing typos
reg_brick_masked[reg_brick_masked == 215] <- 1 #anthropic - fixing typos
reg_brick_masked[reg_brick_masked > 700 & reg_brick_masked < 800] <- NA # remove values that account for misclassification
reg_brick_masked[reg_brick_masked > 400 & reg_brick_masked < 500] <- NA # remove values that account for urban areas/misclassification

# now, we remove pixels that show unflagged moments of regrowth or suppression.
# Regrowth is supposed to be classified as 503, and by looking at the data
# it’s clear that some moments aren’t identified as regrowth correctly by the algorithm,
# looking like ...100—100—100—303—303—303...
#To avoid issues with misclassification, moments that were not designated as proper regrowth by Mapbiomas were removed

# make one empty raster, added to the beginning and the end of the brick.
# SECONDARY IS 0 AND ANTHROPIC IS 1;
r1 <- terra::rast(ncol = ncol(reg_brick_masked), nrow = nrow(reg_brick_masked), ext = ext(reg_brick_masked))
r1 <- terra::setValues(r1, NA)
tmp = c(r1, reg_brick_masked)
tmp2 = c(reg_brick_masked,r1)
tmp3 = tmp-tmp2
# after subtracting one from the other, we get
# -1 is evidence of unflagged suppression
# 1 is unflagged regrowth
tmp3[tmp3 == 1 & tmp3 == -1] <- NA
tmp4 <- tmp3[[2:32]] # remove all-NA columns
tmp5 <- mask(tmp4, sum(tmp4)) # remove all pixels that contain AT LEAST one NA value
# now, all layers of tmp5 should contain only pixels that are complete with data,
# and show no unflagged regrowth/suppression in its history

regrowth_cleaned <- mask(reg_brick_masked, tmp5[[1]])

# calculate age of forest, having as reference year in which regrowth was detected last (2019, for when we have biomass data)
# the which.lyr function returns the FIRST encountering of the value given.
# since in this case we want the most recent regrowth event, we flip the order of the layers in the raster.
regrowth_cleaned_flipped <- regrowth_cleaned [[ c(rev(order(names(regrowth_cleaned)))) ]]
forest_age <- which.lyr(regrowth_cleaned_flipped == 503)

writeRaster(forest_age, paste0('./model_ready_rasters/', location, '_forest_age.tif'))

# 1.72908 % of the pixels are nonNA, meaning they are currently regrowing forests without any classification issues.
