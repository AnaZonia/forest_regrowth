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
library(tidyverse)
library(tidyverse)
library(plyr) # for function ldply
## Brazil shapefile mask
library(maptools)  ## For wrld_simpl
library(pbapply) #progress bar for apply family of functions
data(wrld_simpl)
setwd("/home/aavila/forest_regrowth")
# sourcing functions
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")

##########  Switches ##########
download_worldclim = FALSE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  Temperature/Precipitation ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
outdir <- "./worldclim_dataframes" # Specify your preferred working directory
location <- '0000000000-0000095232'

if (download_worldclim == T){
  vars = c("tmin_", "tmax_", "prec_")
  ranges = c("1980-1989", "1990-1999", "2000-2009", "2010-2018")

  url = paste0("http://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_", do.call(paste0, expand.grid(vars, ranges)), ".zip")
  zip <- file.path(outdir, basename(url))

  for (i in 1:6){
    download.file(url[i], zip[i])
  }

  # get all the zip files
  zipF <- list.files(path = outdir, pattern = "*.zip", full.names = TRUE)
  # unzip all your files
  ldply(.data = zipF, .fun = unzip, exdir = outdir)
  #remove all datapoints before landsat (1985)
  file.remove(list.files(path = outdir, pattern = "1980|1981|1982|1983|1984", full.names = TRUE))
}

mk_list <- function(x)  {as.list(intersect(list.files(path = outdir, pattern = "*.tif", full.names = TRUE),
                                           list.files(path = outdir, pattern = x, full.names = TRUE)))}

vars <- c("prec", "tmax", "tmin")

climate_BRA <- function (var){
  list_clim <- mk_list(var)
  # read the files
  raster_clim <- lapply(list_clim, raster) # this is correct
  # if we are subsetting the data into our region of interest
  BRA <- subset(wrld_simpl, NAME=="Brazil") # get Brazil mask
  rstr <- pbapply::pblapply(raster_clim, terra::crop, BRA) # subsects all rasters to area of interest
  rstr <- pbapply::pblapply(rstr, terra::mask, BRA) # subsects all rasters to area of interest
  rstr <- brick(rstr)
return(rstr)}

prec_BRA <- climate_BRA(vars[1])
tmax_BRA <- climate_BRA(vars[2])
tmin_BRA <- climate_BRA(vars[3])

# processes the dataframes into yearly climate data
climate_yearly = function(tif, fun){
  dates = substr(names(tif), start = 17, stop = 24)   # makes column names only "yyyy.mm"
  dates
  #sum layers, get variance
  indices <- rep(1:34,each=12)
  yearly_data <- stackApply(tif, 1:34, fun = fun)
  names(yearly_data) = c(1985:2018)
  return(yearly_data)
}

temp_BRA = (tmax_BRA + tmin_BRA) / 2

prec_BRA_mean <- climate_yearly(prec_BRA, mean)
prec_BRA_sd <- climate_yearly(prec_BRA, sd)
temp_BRA_mean <- climate_yearly(temp_BRA, mean)
temp_BRA_sd <- climate_yearly(temp_BRA, sd)

writeRaster(prec_BRA_mean, filename='prec_BRA_mean.tif', format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(prec_BRA_sd, filename='prec_BRA_sd.tif', format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(temp_BRA_mean, filename='temp_BRA_mean.tif', format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
writeRaster(temp_BRA_sd, filename='temp_BRA_sd.tif', format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

#resolution is about 4.5 km.

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