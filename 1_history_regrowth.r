####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Extracting regrowth information for the Brazilian Amazon
# Ana Avila - Jan 2023
# Intakes:
#   Raster mask for regrowth history
#   Raw regrowth rasters
# Outputs:
#   regrowth dataframe (full history for all pixels showing regrowth)
#   forest_age dataframe (time of last observed regrowth)
####################################################################

library(sf)
library(terra) # handling spatial data
library(stringr)
#install.packages('pbapply') #progress bar for apply family of functions
library(raster) # because we need to use the getValues function
setwd("/home/aavila/forest_regrowth")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r") # sourcing functions

path <- './mapbiomas/regrowth_rasters'
files <- list.files(path)
locations <- str_sub(files, end = -10)
locations <- unique(locations) #gets all locations currently already processed into filtered, regrowth-only rasters

#for (location in locations){}
# since not all masks are ready and computing this is time consuming, we start with an example location:
location = '0000000000-0000095232'
files_tmp <- list.files(path = './mapbiomas/regrowth_raw', pattern=location, full.names=TRUE)   # obtain paths for all files for that location
files_tmp <- sort(files_tmp)
regrowth_mask <- rast(paste0('./mapbiomas/regrowth_masks/', location, '_regrowth_mask.tif'))

tmp_rasters <- lapply(files_tmp, rast)
reg_brick <- rast(tmp_rasters)

reg_brick_cropped <- crop(reg_brick, ext(regrowth_mask))
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
# select only years that have regrowth that hasn't been suppressed.

regrowth_instances <- which.lyr(reg_brick_masked == 503)
regrowth_last_instance <- where.max(regrowth_instances)
suppression_instances <- which.lyr(600 <= reg_brick_masked & reg_brick_masked < 700)
suppression_last_instance <- where.max(suppression_instances)

# removing all instances of suppression happening after regrowth
delta_regrowth <- regrowth_last_instance - suppression_last_instance
delta_regrowth[delta_regrowth < 0] <- NA
delta_regrowth[delta_regrowth > 0] <- 1
reg_brick_masked <- mask(reg_brick_masked, delta_regrowth)

##################################################################################################

# this removes any instance of -600 coming after -500
reg_brick_masked[reg_brick_masked < 400 & reg_brick_masked >= 300] <- 0 #secondary
reg_brick_masked[100 <= reg_brick_masked & reg_brick_masked < 200] <- 1 #anthropic
reg_brick_masked[reg_brick_masked == 515] <- 1 #anthropic - fixing typos
reg_brick_masked[reg_brick_masked == 215] <- 1 #anthropic - fixing typos
reg_brick_masked[reg_brick_masked > 700 & reg_brick_masked < 800] <- NA # remove values that account for misclassification
reg_brick_masked[reg_brick_masked > 400 & reg_brick_masked < 500] <- NA # remove values that account for urban areas/misclassification

# now, we remove pixels that show unflagged moments of regrowth or repression
r1 <- terra::rast(ncol = ncol(reg_brick_masked), nrow = nrow(reg_brick_masked))
r1 <- terra::setValues(r1, NA)
tmp = c(r1, regrowth_last_instance)
tmp2 = c(regrowth_last_instance,r1)
tmp3 = tmp-tmp2

# -1 is evidence of unflagged repression
# 1 is unflagged regrowth
tmp3[tmp3 == 1] <- NA
tmp3[tmp3 == -1] <- NA
tmp3 = tmp3[,2:(nlyr(tmp)-1)]

# calculate age of forest, having as reference year in which regrowth was detected last (2019, for when we have biomass data)
regrowth_cleaned <- cbind(regrowth_cleaned, 'forest_age' = max(regrowth_cleaned$last_regrowth)-regrowth_cleaned$last_regrowth)

writeRaster(regrowth_cleaned, paste0(location, '_regrowth_history.rds'))
