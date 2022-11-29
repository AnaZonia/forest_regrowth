####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the district of Paragominas, in the Western Amazon.
# Ana Avila - 2022
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

setwd("/home/aavila/Documents/forest_regrowth")


####################################################################
########## FUNCTIONS ##########
####################################################################

# Finds utm zone from longitude - allows to analyze data spanning multiple UTM zones
# This function is necessary for LongLatToUTM (below)
# It intakes:
# Longitude (numeric)
# It outputs:
# UTM zone (numeric)
long2UTMzone = function(long) {
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
LongLatToUTM = function(x,y){
  xy = data.frame(x = x, y = y)
  xy$zone = long2UTMzone(x)

  # split the dataframe by UTM zones
  list_zones = split(xy, xy$zone)

  res_list = list()

  #must convert coodinates separately for different zones
  for (i in 1:length(list_zones)){
    z = list_zones[[i]][1,ncol(list_zones[[i]])] #obtain zone value
    coordinates(list_zones[[i]]) = c("x", "y")
    proj4string(list_zones[[i]]) = CRS("+proj=longlat +datum=WGS84") #convert to spatial object
    # EPSG code calculated for the southern hemisphere as 32700+UTM zone
    # add converted coordinates back into the list
    # obtain list of SpatialObjects, one per UTM zone, with the UTM coordinates.
    res_list = append(res_list, spTransform(list_zones[[i]], CRS(paste("+proj=utm +zone=", z, " +init=epsg:327", z, sep=''))) )
  }

  #convert SpatialObjects into data frames
  res_list = lapply(res_list, as.data.frame)

  #if our data spans more than one UTM zone, res_list will have more than one element.
  #in this case, we unite them all into a single dataframe with zone, x and y columns.
  if (length(res_list) != 1){
    for (i in 2:length(res_list)){
    res_list[[1]] = rbind(res_list[[1]], res_list[[i]])
    }
  }
  
  #convert all coordinates to the nearest multiple of 30 (nearest pixel present in Landsat resolution)
  res_list[[1]]$x = round( res_list[[1]]$x/30 ) * 30
  res_list[[1]]$y = round( res_list[[1]]$y/30 ) * 30

  result = res_list[[1]][ order(as.numeric(row.names(res_list[[1]]))), ]

  #returns dataframe with zone, x and y columns.
  return(result)
}

# Unifies all dataframes in a list into a single dataframe with one column per year.
# Useful when you have a list with one dataframe per moment in time.
# It intakes:
# List of dataframes
# It outputs:
# One combined dataframe containing the THIRD column of every dataframe in the list
df_merge = function(df){
  for (i in 2:length(df)){
    df[[1]] = cbind(df[[1]], df[[i]][3])
  }
  return(df[[1]])
}

# reduce date from "%Y-%m-%d %H:%M:%S" format into just the year
extract_year = function(df, in_format, out_format){
  df$date = as.POSIXct(df$date, format = in_format)
  df$date = format(df$date, format=out_format) 
  return(df)
}

# find last instance of a value in a dataframe, and returns the year (column name) of that occurrence
find_last_instance = function(df, fun){
  instances = sapply(apply(df[3:(ncol(df)-3)], 1, fun),names) # returns rowname/index as well as years flagged as FIRE events
  last = lapply(lapply(instances, unlist), max) # select for only the most recent observed regrowth
  last_df = as.data.frame(unlist(last))
  last_df = as.data.frame(lapply(last_df,as.numeric))
  last_df[is.na(last_df)] = 0
  return(last_df)
}

####################################################################
########## SWITCHES ##########
####################################################################
import_mapbiomas = T
import_clim == T
import_dubayah = F
import_santoro = F
import_potapov = F

####################################################################
########## EXTRACTING DATA ##########
####################################################################

##########  REGROWTH/DEFORESTATION, FIRE AND LAND USE ##########


#specifying coordinates of interest - where would you like to crop the dataframe to?

# EXTRACTING DATA FROM MULTIPLE REGIONS OF THE COUNTRY (larger scale)

files <- list.files(path = './mapbiomas/regrowth_amazon')
locations <- str_sub(files, start= -25, end = -5)
locations <- unique(locations)
dir.create('regrowth_dataframes')

# the Amazon is divided in 12 parts, each with their own identifier
# each location is an identifier.
for (location in locations){
  files_tmp <- list.files(path = './mapbiomas/regrowth_amazon', pattern=location, full.names=TRUE)   # obtain paths for all files for that location

  dir.create(location)

  # obtain the raster file for all years within that location
  regrowth_list = c()
  for (i in 1:length(files)){
    regrowth_list[[i]] = raster(files[i])
    print(i)
  }

  # NOTE: MAY NEED TO ORDER THE FILE NAMES BY YEAR

  # to make processing lighter, subset only the pixels that have shown regrowth history.
  # here, we are (1) making a mask registering all regrowth moments and (2) subsetting rasters based on that mask.
  # a regrowth moment is flagged with the value "503", therefore:
  for (i in 1:length(regrowth_list)){
    regrowth_list[[i]][regrowth_list[[i]]!=503] <- NA # only leave behind values not including the regrowth moment
    writeRaster(tmp_dfs[[i]], file.path(location, paste0(c(1984+i), ".tif"))) # save rasters with the year on the folder created for each location.
  }

  # the reason these files are being stored in the machine is to avoid losing all progress in case of a disconnection of the ssh,
  # as well as avoid having to redo everything in case there is an error that could be fixed later.

  # create a raster stack with one layer per year, for this location. This will be a large stack with only 0 or 503.
  stacked_years = raster(paste0(location, '1985.tif'))
  for (i in 1986:2018){
    tmp = raster(paste0(i, '.tif'))
    tst = stack(stacked_years, tmp)
  }

  #stacked files are merged to become one 
  regrowth_mask = merge(stacked_years)

  masked = lapply(regrowth_list, mask, regrowth_mask)

  stacked_history = stack(masked)

  # convert into dataframe for further manipulation
  df_tst = as.data.frame(stacked_history, xy = T, na.rm = T)
  saveRDS(df_tst, paste0(, ""))

  #replace function from dplyr should work here
  colnames(df_tst) = c(sub('mapbiomas.brazil.collection.60.', "", colnames(df_tst[,1:c(ncol(df_tst)-2)])), "lon", "lat")
  df_tst <- df_tst[ , order(names(df_tst))]
  #missing 1990, 1994, 1999



}



saveRDS(df_tst, "df_tst.rds")

writeRaster(tst, "masked_merged.tif")

tst = raster("masked_merged.tif")

writeRaster(tst_mrg, "merged_years.tif")



##########  TEMPERATURE AND RAINFALL ##########

if (import_clim == T){

  climate_data_import = function(df){
    colnames(df) = str_sub(colnames(df), start= -7)   # makes column names only "yyyy.mm"
    result = t(apply(df[3:ncol(df)], 1, tapply, gl(34, 12), mean))
    colnames(result) = c(1985:2018)
    return(as.data.frame(cbind(df[,1:2], result)))
  }

  df_prec = readRDS('df_prec.rds')
  df_tmin = readRDS('df_tmin.rds')
  df_tmax = readRDS('df_tmax.rds')

  prec = climate_data_import(df_prec)
  tmax = climate_data_import(df_tmax)
  tmin = climate_data_import(df_tmin)
  temp = (tmax + tmin) / 2

  climate_data_cleanup = function(df){
    df$mean = colMeans(df[3:ncol(df)])
    df = cbind(df, LongLatToUTM(df$x, df$y))
    df = df[,3:ncol(df)]
    df$xy = paste0(df$zone, df$x, df$y)
    return(df)
  }

  prec = climate_data_cleanup(prec)
  temp = climate_data_cleanup(temp)

  #resolution is about 4.5 km.
  #Need to add temperature data for ALL cells within 4.5km of center point.

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

  saveRDS(df_prec, "df_prec.rds")
  saveRDS(df_tmax, "df_tmax.rds")
  saveRDS(df_tmin, "df_tmin.rds")
}

##########  SOIL TYPE ##########

soil <- sf::st_read( 
  dsn= paste0(getwd(),"/soil/") , 
  layer="DSMW"
)

brazil_soil = soil[soil$COUNTRY == "BRAZIL",]

# SNUM - Soil mapping unit number
# 

##########  BIOMASS ##########


#The CRAN version:
install.packages("rGEDI")

#The development version:
library(devtools)
devtools::install_github("carlos-alberto-silva/rGEDI", dependencies = TRUE)

# loading rGEDI package
library(rGEDI)



# Dubayah et al 2022 -> GEDI L4A Footprint Level Aboveground Biomass Density (Mg/ha)
# 1km resolution, 2019-2021
if (import_dubayah == T){
  biomass = readRDS("./biomass_dubayah.rds")
}else{
  
  lat_lon_oi = c(ymax, ymin, xmin, xmax)  #lat and lon of interest

  GEDI_download = l4_download(
    lat_lon_oi,
    outdir = "./GEDI",
    from = "2020-01-01",
    to = "2020-07-31",
    just_path = F
  )

  filepaths = paste0('./GEDI/', list.files("./GEDI"))
  GEDI_list = lapply(filepaths, l4_get, just_colnames = F)
  GEDI_list = lapply(GEDI_list, subset, select=c("date", "lat_lowestmode", "lon_lowestmode", "agbd_se", "agbd"))
  GEDI_list = lapply(GEDI_list, extract_year, in_format = "%Y-%m-%d %H:%M:%S", out_format = "%Y")

  select_range = function(df){
    df = df[ymin < lat_lowestmode & lat_lowestmode < ymax & xmin < lon_lowestmode & lon_lowestmode < xmax,]
    return(df)
  }

  GEDI_list = lapply(GEDI_list, select_range)

  for (i in 2:length(GEDI_list)){
  GEDI_list[[1]] = rbind(GEDI_list[[1]], GEDI_list[[i]])
  }

  biomass = GEDI_list[[1]]

  biomass = cbind(biomass, LongLatToUTM(biomass$lon_lowestmode, biomass$lat_lowestmode))

  saveRDS(biomass, "biomass_dubayah.rds")
}

# Santoro et al 2018 data -> GlobBiomass ESA (Mg/ha)
# 100m resolution, 2010
if (import_santoro == T){
  biomass = readRDS("biomass_santoro_Brazil.rds")
}else{
  biomass1 = raster("N00W060_agb.tif")
  biomass2 = raster("N00W100_agb.tif")
  biomass3 = raster("N40W060_agb.tif")

  biomass = merge(biomass1, biomass2, biomass3)

  ## crop and mask
  r2 <- crop(biomass, extent(BRA))
  r3 <- mask(r2, BRA)

  biomass = as.data.frame(biomass, xy = TRUE)

  #biomass = biomass[ymin < biomass$y & biomass$y < ymax & xmin < biomass$x & biomass$x < xmax,]

  biomass = cbind(biomass, LongLatToUTM(biomass$x, biomass$y))

  colnames(biomass) = c('lon', 'lat', 'agbd', 'zone', 'x', 'y')
  saveRDS(biomass, "biomass_santoro.rds")

}

# Potapov et al 2020 -> GLAD Forest Canopy Height (m)
# 30m resolution, 2019
if (import_potapov == T){
  biomass = readRDS("Forest_height_2019_Brazil.rds")
}else{
  biomass = raster("Forest_height_2019_SAM.tif")

  ## crop and mask
  r2 <- crop(biomass, extent(BRA))
  r3 <- mask(r2, BRA)

  #writeRaster(r3, "Forest_height_2019_Brazil.tif")

}