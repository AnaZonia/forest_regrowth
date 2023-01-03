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
library(raster) #  handling spatial data
library(terra) # handling spatial data
library(geodata) # to extract worldclim with getData
library(sp) # to extract worldclim with getData
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

regrowth_list <- c()
for (i in 1:length(files_tmp)){
  regrowth_list <- c(regrowth_list, raster(files_tmp[i]))
  print(i)
}

# now, we use the regrowth mask to make a regrowth stack with all history per location,
#containing data only on pixels that show any regrowth event in their history.
regrowth_mask <- raster(paste0('./mapbiomas/regrowth_masks/', location, '_mask.tif'))
masked <- pbapply::pblapply(regrowth_list, terra::mask, regrowth_mask) # mask all raw files, leaving behind only the cells containing any regrowth instance.
stacked_history <- stack(masked) # now, we have a stack with all regrowth/deforestation history for a location.
#writeRaster(stacked_history, "0000000000-0000095232_regrowth.tif")
# convert into dataframe for further manipulation
# df_history <- as.data.frame(stacked_history, xy = T, na.rm = T)

# This stack gets quite large. as.data.frame and getValues both break when handling it because of its size.
# even writeRaster takes 30min+ to work with a stack this size.
# My solution was to work with a smaller part at a time.
e = extent(-48.31863, -43.99998, -3.2823093, -0.5377764)
stacked_history1 <- terra::crop(stacked_history, e)

convert_history <- getValues(stacked_history1)
convert_history <- data.frame(cell = 1:length(convert_history), value = convert_history)
convert_history <- na.omit(convert_history)
convert_history[,c("x","y")] <- xyFromCell(stacked_history1, convert_history$cell)

# cleaning column names (which come out long and messy)
# we want each column name to contain only the year of the data.
location_colname <- paste0('.', gsub('-', '.', location))
subs <- c('mapbiomas.brazil.collection.60.', location_colname)
for (cut in subs){
  colnames(convert_history) <- c(sub(cut, "", colnames(convert_history[,1:c(ncol(convert_history)-2)])), "lon", "lat")
}

# adding UTM coordinates in case it's needed
saveRDS(convert_history, file.path(paste0('./mapbiomas/dataframes/', location, '_regrowth.rds')))
#regrowth = readRDS('0000000000-0000095232_regrowth.rds')

###################################
########## DATA CLEANING ##########
###################################

setwd("/home/aavila/forest_regrowth/mapbiomas/dataframes")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r") # sourcing functions
regrowth = readRDS('0000000000-0000095232_regrowth_2.rds')
head(regrowth)

# INDEX #
# 100 anthropic
# 200 mature
# 300 secondary
# 400 deforestation
# 500 regrowth
# 600 secondary suppression

#if (import_santoro == T){    # since santoro data is available for 2010 rather than 2019,
#  regrowth = cbind(regrowth[,1:25], regrowth[,(ncol(regrowth)-2):ncol(regrowth)])
#  fire = cbind(fire[,1:28], fire[,(ncol(fire)-2):ncol(fire)])
#  lulc = cbind(lulc[,1:28], lulc[,(ncol(lulc)-2):ncol(lulc)])
#}

# select only years that have regrowth that hasn't been suppressed.
regrowth_last_instance <- find_last_instance(regrowth, function(x) which(x == 503))
colnames(regrowth_last_instance) <- "last_regrowth"

suppression_last_instance = find_last_instance(regrowth, function(x) which(600 <= x & x < 700))
colnames(suppression_last_instance) = "last_suppression"

# find instances in which the last regrowth detected hasn't been followed by suppression
regrowth_unsuppressed <- cbind(regrowth[(ncol(regrowth)-1):ncol(regrowth)],regrowth_last_instance)
regrowth_unsuppressed <- subset(regrowth_unsuppressed, regrowth_unsuppressed$last_regrowth-suppression_last_instance > 0)
regrowth_unsuppressed <- cbind(regrowth_unsuppressed[,(ncol(regrowth_unsuppressed)-2):ncol(regrowth_unsuppressed)], 'forest_age' = max(regrowth_unsuppressed$last_regrowth)-regrowth_unsuppressed$last_regrowth)

##################################################################################################

# this removes any instance of -600 coming after -500
regrowth2 = regrowth[rownames(regrowth) %in% rownames(regrowth_unsuppressed), ] 
regrowth2[regrowth2 < 400 & regrowth2 >= 300] = 0 #secondary
regrowth2[100 <= regrowth2 & regrowth2 < 200] = 1 #anthropic
regrowth2[regrowth2 == 515] = 1 #anthropic - fixing typos
regrowth2[regrowth2 == 215] = 1 #anthropic - fixing typos

regrowth2[regrowth2 > 700 & regrowth2 < 800] <- NA # remove values that account for misclassification
regrowth2[regrowth2 > 400 & regrowth2 < 500] <- NA # remove values that account for urban areas/misclassification
regrowth2 = regrowth2[complete.cases(regrowth2),]

# now, we remove pixels that show unflagged moments of regrowth or repression. (shoutout to Jacqueline Oehri for method!)
tmp = cbind(NA, regrowth2[,3:(ncol(regrowth)-3)])
tmp2 = cbind(regrowth2[,3:(ncol(regrowth)-3)],NA)
tmp3 = tmp-tmp2

# -1 is evidence of unflagged repression
# 1 is unflagged regrowth
tmp3[tmp3 == 1] <- NA
tmp3[tmp3 == -1] <- NA
tmp3 = tmp3[,2:(ncol(tmp)-1)]
tmp3 = tmp3[complete.cases(tmp3),]

#selects for cleaned rows
regrowth <- regrowth[rownames(regrowth) %in% rownames(tmp3), ]
regrowth_last_instance <- find_last_instance(regrowth, function(x) which(x == 503))
head(regrowth_last_instance)
colnames(regrowth_last_instance) <- "last_regrowth"

regrowth_cleaned <- cbind(regrowth[(ncol(regrowth)-1):ncol(regrowth)],regrowth_last_instance) # create dataframe with years and last instance of growth
# calculate age of forest, having as reference year in which regrowth was detected last (2019, for when we have biomass data)
regrowth_cleaned <- cbind(regrowth_cleaned, 'forest_age' = max(regrowth_cleaned$last_regrowth)-regrowth_cleaned$last_regrowth)
regrowth_cleaned <- cbind(regrowth_cleaned, LongLatToUTM(regrowth_cleaned$lon, regrowth_cleaned$lat)) # add UTM coordinates as new columns
regrowth_cleaned$xy <- paste0(regrowth_cleaned$zone, regrowth_cleaned$x, regrowth_cleaned$y) # create unique identifier

saveRDS(regrowth_cleaned, paste0(location, '_forest_age.rds'))
