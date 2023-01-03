####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# Extracting biomass and regrowth information for the entire Amazon.
# Ana Avila - Dec 2022
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

setwd("/home/aavila/Documents/forest_regrowth")


####################################################################
##########              SWITCHES             #######################
####################################################################
regrowth <- FALSE 
lulc <- FALSE 
fire <- FALSE 
climate


####################################################################
##########              FUNCTIONS            #######################
####################################################################

# Finds utm zone from longitude - allows to analyze data spanning multiple UTM zones
# This function is necessary for LongLatToUTM (below)
# It intakes:
# Longitude (numeric)
# It outputs:
# UTM zone (numeric)
long2UTMzone <- function(long) {
  ## Function to get the UTM zone for a given longitude
  (floor((long + 180)/6) %% 60) + 1
}

# Converts coodinates from lat/long to UTM.
# this allows merging dataframes with different coordinates.
# It converts all points to to 30m resolution by finding the nearest coordinate multiple of 30.
# It intakes:
# x and y = longitude and latitude, respectively (numeric or vector)
# It outputs:
# Dataframe with zone, x and y columns in UTM format.
LongLatToUTM <- function(x,y){
  xy <- data.frame(x = x, y = y)
  xy$zone <- long2UTMzone(x)

  # split the dataframe by UTM zones
  list_zones <- split(xy, xy$zone)

  res_list <- list()

  #must convert coodinates separately for different zones
  for (i in 1:length(list_zones)){
    z <- list_zones[[i]][1,ncol(list_zones[[i]])] #obtain zone value
    coordinates(list_zones[[i]]) <- c("x", "y")
    proj4string(list_zones[[i]]) <- CRS("+proj=longlat +datum=WGS84") #convert to spatial object
    # EPSG code calculated for the southern hemisphere as 32700+UTM zone
    # add converted coordinates back into the list
    # obtain list of SpatialObjects, one per UTM zone, with the UTM coordinates.
    res_list <- append(res_list, spTransform(list_zones[[i]], CRS(paste("+proj=utm +zone=", z, " +init=epsg:327", z, sep=''))) )
  }

  #convert SpatialObjects into data frames
  res_list <- lapply(res_list, as.data.frame)

  #if our data spans more than one UTM zone, res_list will have more than one element.
  #in this case, we unite them all into a single dataframe with zone, x and y columns.
  if (length(res_list) != 1){
    for (i in 2:length(res_list)){
    res_list[[1]] <- rbind(res_list[[1]], res_list[[i]])
    }
  }
  
  #convert all coordinates to the nearest multiple of 30 (nearest pixel present in Landsat resolution)
  res_list[[1]]$x <- round( res_list[[1]]$x/30 ) * 30
  res_list[[1]]$y <- round( res_list[[1]]$y/30 ) * 30

  result <- res_list[[1]][ order(as.numeric(row.names(res_list[[1]]))), ]

  #returns dataframe with zone, x and y columns.
  return(result)
}

# Unifies all dataframes in a list into a single dataframe with one column per year.
# Useful when you have a list with one dataframe per moment in time.
# It intakes:
# List of dataframes
# It outputs:
# One combined dataframe containing the THIRD column of every dataframe in the list
df_merge <- function(df){
  for (i in 2:length(df)){
    df[[1]] <- cbind(df[[1]], df[[i]][3])
  }
  return(df[[1]])
}

# reduce date from "%Y-%m-%d %H:%M:%S" format into just the year
extract_year <- function(df, in_format, out_format){
  df$date <- as.POSIXct(df$date, format = in_format)
  df$date <- format(df$date, format = out_format) 
  return(df)
}

# find last instance of a value in a dataframe, and returns the year (column name) of that occurrence
find_last_instance <- function(df, fun){
  instances <- sapply(apply(df[3:(ncol(df)-3)], 1, fun),names) # returns rowname/index as well as years flagged as FIRE events
  last <- lapply(lapply(instances, unlist), max) # select for only the most recent observed regrowth
  last_df <- as.data.frame(unlist(last))
  last_df <- as.data.frame(lapply(last_df,as.numeric))
  last_df[is.na(last_df)] <- 0
  return(last_df)
}

# makes dataframe from large raster files, substituting as.data.frame()
# assumes cells in the raster start from 1 for correct assignment of coordinates.
df_from_raster <- function(raster){
  bm_test <- getValues(raster)
  bm_test <- data.frame(cell = 1:length(bm_test), value = bm_test)
  bm_test <- na.omit(bm_test)
  bm_test[,c("x","y")] <- xyFromCell(raster, bm_test$cell)
  return(bm_test)
}

# imports raw regrowth data and returns the rasters in a list.
import_regrowth <- function(location){
  files_tmp <- list.files(path = './mapbiomas/regrowth_raw', pattern=location, full.names=TRUE)   # obtain paths for all files for that location
  files_tmp <- sort(files_tmp)

  regrowth_list <- c()
  for (i in 1:length(files_tmp)){
    regrowth_list <- c(regrowth_list, raster(files_tmp[i]))
    print(i)
  }

  return(regrowth_list)
}

####################################################################
##########               BODY                #######################
####################################################################



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  REGROWTH/DEFORESTATION ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#specifying coordinates of interest - where would you like to crop the dataframe to?

# EXTRACTING DATA FROM MULTIPLE REGIONS OF THE COUNTRY (larger scale)
path <- './mapbiomas/regrowth_raw'
files <- list.files(path)
#newname <- sub('none-','', files) ## making names easier to read, standardizing the names
#file.rename(file.path(path,files), file.path(path, newname)) ## renaming it.

locations <- str_sub(files, start= -25, end = -5)
locations <- unique(locations)
#locations <- locations[1:length(locations)]


# the Amazon is divided in 12 parts, each with their own identifier (a location number).
# here we go location by location:
  # import raw data
  # make a mask with only areas that show regrowth - this mask will be used to reduce the size of files we're handling
  # save these rasters in the ./regrowth_rasters directory
for (i in 1:length(locations)){ 
  location <- '0000031744-0000063488'
  regrowth_list <- import_regrowth(location)

  # obtain the raster file for all years within that location.
  # to make processing lighter, subset only the pixels that have shown regrowth history.
  # here, we are (1) making a mask registering all regrowth moments and (2) subsetting rasters based on that mask.
  # a regrowth moment is flagged with the value "503", therefore:
  for (i in 1:1){
    print(i)
    print(Sys.time())
    regrowth_list[[i]][regrowth_list[[i]]!=503] <- NA # only leave behind values not including the regrowth moment
    writeRaster(regrowth_list[[i]], file.path(paste0('./mapbiomas/regrowth_rasters/', location, '_', c(1987+i), ".tif"))) # save rasters with the year on the folder created for each location.
  }

}

# the reason these files are being stored in the machine is to avoid losing all progress in case of a disconnection of the ssh,
# as well as avoid having to redo everything in case there is an error that could be fixed later by editing the files.
# create a raster stack with one layer per year, for this location.
# this stack will be merged to make the regrowth-only mask.
path <- './mapbiomas/regrowth_rasters'
files <- list.files(path)
locations <- str_sub(files, end = -10)
locations <- unique(locations) #gets all locations currently already processed into filtered, regrowth-only rasters

#for (location in locations){
  location = "0000000000-0000095232"
  filepaths <- paste0("/home/aavila/Documents/forest_regrowth/mapbiomas/regrowth_rasters/", location, '_', c(1988:2019), '.tif')
  mask_raster_list <- lapply(filepaths, raster)
  mask_raster_stack <- stack(mask_raster_list)
  regrowth_mask <- merge(mask_raster_stack)
  regrowth_mask <- writeRaster(regrowth_mask, file.path(paste0('./mapbiomas/regrowth_masks/', location, '_mask.tif'))) # mask made and saved.

  # now, we use the regrowth mask to make a regrowth stack with all history per location,
  #containing data only on pixels that show any regrowth event in their history.
  #regrowth_mask <- raster('./mapbiomas/regrowth_masks/0000000000-0000095232_mask.tif')
  regrowth_list <- import_regrowth(location)
  masked <- pbapply::pblapply(regrowth_list, terra::mask, regrowth_mask)
  stacked_history <- stack(masked) # we have a stack with all history for a location.

  # This stack gets quite large. as.data.frame and getValues both break when handling it because of its size.
  # even writeRaster takes 30min+ to work with a stack this size.
  # My solution was to break it in half and work with smaller parts at a time.
  
  #writeRaster(stacked_history, "0000000000-0000095232_regrowth.tif")

  # lulc = readRDS('./test_files/lulc.rds')
  # lulc2 = lulc[,c(1:(ncol(lulc)-3))]
  # x <- c(min(lulc2$lon), max(lulc2$lon))
  # y <- c(min(lulc2$lat), max(lulc2$lat))
  # e = extent(x, y)
# > range(lulc$lat)
# [1] -3.2823093 -0.5377764
# > range(lulc$lon)
# [1] -48.32644 -43.99998
#  xmin, xmax, ymin, ymax

  stacked_history1 <- terra::crop(stacked_history, e)

  
  convert_history <- getValues(stacked_history1)
  convert_history <- data.frame(cell = 1:length(convert_history), value = convert_history)
  convert_history <- na.omit(convert_history)
  convert_history[,c("x","y")] <- xyFromCell(stacked_history1, convert_history$cell)
  saveRDS(convert_history, file.path(paste0('./mapbiomas/dataframes/', location, '_regrowth.rds')))


  stacked_history2 <- terra::crop(stacked_history, e)


  # convert into dataframe for further manipulation
  # df_history <- as.data.frame(stacked_history, xy = T, na.rm = T)
  # cleaning column names
  location_colname <- paste0('.', gsub('-', '.', location))
  subs <- c('mapbiomas.brazil.collection.60.', location_colname)
  for (cut in subs){
    colnames(df_tst) <- c(sub(cut, "", colnames(df_tst[,1:c(ncol(df_tst)-2)])), "lon", "lat")
  }
  df_history <- cbind(df_history, LongLatToUTM(df_history$lon, df_history$lat))

#}


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


location = locations[4]
# the Amazon is divided in 12 parts, each with their own identifier
# each location is an identifier.
#for (location in locations){
files_tmp <- list.files(path = './mapbiomas/fire_amazon', pattern=location, full.names=TRUE)   # obtain paths for all files for that location
files_tmp <- sort(files_tmp)

#dir.create(paste0('./regrowth_dataframes/', location))

fire_list <- c()
for (i in 4:length(files_tmp)){ # since regrowth info is only available 1988 onwards
  fire_list <- c(fire_list, raster(files_tmp[i]))
  print(i)
}

  # obtain the raster file for all years within that location

# to make processing lighter, subset only the pixels that have shown regrowth history.



fire_mask <- raster('regrowth_mask.tif')

fire_list <- lapply(fire_list, crop, fire_mask)

fire_mask <- crop(fire_mask, fire_list[[1]])

masked <- lapply(fire_list, mask, fire_mask)

fire_history <- stack(masked)

df_tst <- as.data.frame(fire_history, xy = T, na.rm = T)








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  Land use/Land cover ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


files <- list.files(path = './mapbiomas/lulc', pattern='\\.tif$', full.names=TRUE)   # obtain paths for all files 
tmp_rasters <- lapply(files, raster)

regrowth_mask <- raster('./mapbiomas/regrowth_masks/0000000000-0000095232_mask.tif')

# if we are subsetting the data into our region of interest
# coord_oi <- c(-48.31876, -43.99984, -3.282444, 5.269697)  #(xmin, xmax, ymin, ymax)
# e <- as(extent(coord_oi), 'SpatialPolygons')
# crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"

tmp_rasters <- pbapply::pblapply(tmp_rasters, terra::crop, extent(regrowth_mask)) # subsects all rasters to area of interest
tmp_rasters2 <-  pbapply::pblapply(tmp_rasters, terra::mask, regrowth_mask) # subsects all rasters to area of interest

stacked_tst <- stack(tmp_rasters2)

writeRaster(stacked_tst, "0000000000-0000095232_lulc.tif")


convert_history <- getValues(stacked_history1)
convert_history <- data.frame(cell = 1:length(convert_history), value = convert_history)
convert_history <- na.omit(convert_history)
convert_history[,c("x","y")] <- xyFromCell(stacked_history1, convert_history$cell)
saveRDS(convert_history, file.path(paste0('./mapbiomas/dataframes/', location, '_regrowth.rds')))





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  Temperature/Precipitation ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (import_clim == F){
  temp <- readRDS('./worldclim_brazil/dataframes/temp.rds')
  prec <- readRDS('./worldclim_brazil/dataframes/prec.rds')
}else{

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


}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  Soil ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


soil <- sf::st_read( 
  dsn= paste0(getwd(),"/soil/") , 
  layer="DSMW"
)

brazil_soil = soil[soil$COUNTRY == "BRAZIL",]

test_coords <- lapply(brazil_soil$geometry, st_coordinates)
test_coords <- lapply(test_coords, as.data.frame)
subset2 <- function(df){return(df[-c(3,4)])}
test_coords <- lapply(test_coords, subset2)

test_coords2 <- lapply(test_coords, cbind())


tst <- lapply(soil$test_coords, as.data.frame)

# SNUM - Soil mapping unit number

