####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Prepares fire history data
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Intakes:
#   Raster mask for regrowth history
#   Raw fire rasters
# Outputs:
#   fire dataframe (full history for all pixels showing yes/no burning detected)
#   burn_history dataframe (time since last observed fire; total number of fires)
####################################################################

library(sf)
library(raster) #  handling spatial data
library(terra) # handling spatial data
library(geodata) # to extract worldclim with getData
library(sp) # to extract worldclim with getData
#library(stringr)
library(tidyverse)
library(plyr)
library(foreach) # for splitting heavy processing (masking, converting)
library(doParallel) # for splitting heavy processing (masking, converting)
## Brazil shapefile mask
library(maptools)  ## For wrld_simpl
library(​data.table​) #for faster reading of csv files with function fread
library(pbapply) #progress bar for apply family of functions
data(wrld_simpl)
BRA <- subset(wrld_simpl, NAME=="Brazil")

setwd("/home/aavila/forest_regrowth/dataframes")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  FIRE ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# EXTRACTING DATA FROM MULTIPLE REGIONS OF THE COUNTRY (larger scale)
path <- './mapbiomas/fire_amazon'
files <- list.files(path)
#newname <- sub('none-','', files) ## making names easier to read, standardizing the names
#file.rename(file.path(path,files), file.path(path, newname)) ## renaming it.
locations <- str_sub(files, start= -25, end = -5)
locations <- unique(locations)
#for (location in locations){}
location = '0000000000-0000095232'

files_tmp <- list.files(path = './mapbiomas/fire_amazon', pattern=location, full.names=TRUE)   # obtain paths for all files for that location
files_tmp <- sort(files_tmp)

fire_list <- c()
for (i in 4:length(files_tmp)){ # since regrowth info is only available 1988 onwards
  fire_list <- c(fire_list, raster(files_tmp[i]))
  print(i)
}

# now, we use the regrowth mask to make a regrowth stack with all history per location,
#containing data only on pixels that show any regrowth event in their history.
regrowth_mask <- raster(paste0('./mapbiomas/regrowth_masks/', location, '_mask.tif'))
# correcting for small differences between the sizes
fire_list <- lapply(fire_list, crop, regrowth_mask)
fire_mask <- crop(fire_mask, fire_list[[1]])

masked <- pbapply::pblapply(fire_list, terra::mask, regrowth_mask) # mask all raw files, leaving behind only the cells containing any regrowth instance.
stacked_history <- stack(masked) # now, we have a stack with all regrowth/deforestation history for a location.
fire <- as.data.frame(fire_history, xy = T, na.rm = T)

location_colname <- paste0('.', gsub('-', '.', location))
subs <- c('mapbiomas.brazil.collection.10.amazonia.', location_colname)
for (cut in subs){
  colnames(fire) <- c(sub(cut, "", colnames(fire[,1:c(ncol(fire)-2)])), "lon", "lat")
}

saveRDS(fire, file.path(paste0('./mapbiomas/dataframes/', location, '_fire.rds')))

fire <- readRDS('0000000000-0000095232_fire.rds')
head(fire)

#count number of total fires
fire$num_fires <- rowSums(fire[3:(ncol(fire)-3)])
fire_last_instance <- find_last_instance(fire, function(x) which(x == 1))
colnames(fire_last_instance) <- "last_burn"
fire <- cbind(fire,fire_last_instance)
fire$last_burn <- max(fire$last_burn) - fire$last_burn
fire$last_burn[fire$last_burn == max(fire$last_burn)] = NA

fire <- cbind(fire, LongLatToUTM(fire$lon, fire$lat)) # add UTM coordinates as new columns
fire$xy <- paste0(fire$zone, fire$x, fire$y)

saveRDS(fire[,c(c(ncol(fire)-7): ncol(fire))], paste0(location, '_burn_history.rds'))
#saveRDS(fire, file.path(paste0('./mapbiomas/dataframes/', location, '_burn_history.rds')))
