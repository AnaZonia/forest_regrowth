####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepares dataframes for environmental variables, country-wide (BRA):
  # -> Precipitation and temperature from raw worldclim data
  # -> Soil types from FAO World Soil Map
# Generates:
  # temp <- readRDS('./worldclim_brazil/dataframes/temp.rds')
  # prec <- readRDS('./worldclim_brazil/dataframes/prec.rds')
  # soil <- readRDS('./soil/soil.rds')
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
# sourcing functions
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  Temperature/Precipitation ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  outdir <- "./worldclim_brazil" # Specify your preferred working directory
  vars = c("tmin_", "tmax_", "prec_")
  ranges = c("1980-1989", "1990-1999", "2000-2009", "2010-2018")

  url = paste0("http://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_", do.call(paste0, expand.grid(vars, ranges)), ".zip")
  zip <- file.path(outdir, basename(url))

  for (i in 1:length(url)){
    download.file(url[i], zip[i])
  }

  # get all the zip files
  zipF <- list.files(path = outdir, pattern = "*.zip", full.names = TRUE)
  # unzip all your files
  ldply(.data = zipF, .fun = unzip, exdir = outdir)
  # get the csv files
  #remove all datapoints before landsat (1985)
  file.remove(list.files(path = outdir, pattern = "1980|1981|1982|1983|1984", full.names = TRUE))

  clim_files <- as.list(list.files(path = outdir, pattern = "*.tif", full.names = TRUE))
  # read the csv files
  raster_clim <- lapply(clim_files, raster)
  
  # if we are subsetting the data into our region of interest
  coord_oi = c(xmin, xmax, ymin, ymax) #this specifies the coordinates of Paragominas.
  e = as(extent(coord_oi), 'SpatialPolygons')
  crs(e) = "+proj=longlat +datum=WGS84 +no_defs"
  df_clim = lapply(raster_clim, crop, e) # subsects all rasters to area of interest

  df_clim = lapply(df_clim, as.data.frame, xy=T)

  #we have 409 entries for prec, tmax and tmin.
  df_prec = df_merge(df_clim[1:408])
  df_tmax = df_merge(df_clim[410:817])
  df_tmin = df_merge(df_clim[819:c(length(df_clim))])

  # saveRDS(df_prec, "df_prec.rds")
  # saveRDS(df_tmax, "df_tmax.rds")
  # saveRDS(df_tmin, "df_tmin.rds")

    # processes the dataframes into yearly climate data
  climate_data_import = function(df){
    colnames(df) = str_sub(colnames(df), start= -7)   # makes column names only "yyyy.mm"
    result = t(apply(df[3:ncol(df)], 1, tapply, gl(34, 12), mean)) #mean annual - data is originally monthly
    colnames(result) = c(1985:2018)
    return(as.data.frame(cbind(df[,1:2], result)))
  }

  # df_prec = readRDS('./worldclim_brazil/dataframes/df_prec.rds')
  # df_tmin = readRDS('./worldclim_brazil/dataframes/df_tmin.rds')
  # df_tmax = readRDS('./worldclim_brazil/dataframes/df_tmax.rds')

  prec = climate_data_import(df_prec)
  tmax = climate_data_import(df_tmax)
  tmin = climate_data_import(df_tmin)
  temp = (tmax + tmin) / 2

  # gives proper names to the dataframes and adds UTM coordinates
  climate_data_cleanup = function(df){
    df$mean = colMeans(df[3:ncol(df)])
    names(df)[names(df) == 'x'] <- 'lon'
    names(df)[names(df) == 'y'] <- 'lat'
    df = cbind(df, LongLatToUTM(df$lon, df$lat))
    df$xy = paste0(df$zone, df$x, df$y)
    return(df)
  }

  saveRDS(climate_data_cleanup(prec), file.path(paste0('./worldclim_brazil/dataframes/prec.rds')))
  saveRDS(climate_data_cleanup(temp), file.path(paste0('./worldclim_brazil/dataframes/temp.rds')))
  #resolution is about 4.5 km.
  #Need to add temperature data for ALL cells within 4.5km of center point.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  Soil ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SNUM - Soil mapping unit number

soil <- sf::st_read(dsn= paste0(getwd(),"/soil/") , layer="DSMW")

brazil_soil = soil[soil$COUNTRY == "BRAZIL",]

test_coords <- lapply(brazil_soil$geometry, st_coordinates) #extract lat and lon
test_coords <- lapply(test_coords, as.data.frame) #create a list of dataframes, one per polygon
subset2 <- function(df){return(df[-c(3,4)])} # clean up results of st_coordinates function
test_coords <- lapply(test_coords, subset2)

# add soil types (in DOMSOI category) to coordinate values
test_coords2 <- cbind(test_coords[[1]],type = brazil_soil$DOMSOI[1]) 
colnames(test_coords2) <- c('lon', 'lat', 'type')
for (i in 2:length(brazil_soil$DOMSOI)){
  print(i)
  df = cbind(test_coords[[i]],type = brazil_soil$DOMSOI[i]) # add the 
  #print(df[1,])
  test_coords2 <- rbind(test_coords2, df)
}

saveRDS(test_coords2, 'soil.rds')
