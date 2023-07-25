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
outdir <- "./worldclim" # Specify your preferred working directory
location <- '0000000000-0000095232'
regrowth <- rast("./model_ready_rasters/0000000000-0000095232_forest_age.tif")

if (download_worldclim == T){
  vars = c("tmin_", "tmax_", "prec_")
  ranges = c("1980-1989", "1990-1999", "2000-2009", "2010-2018")

  url = paste0("http://biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_", do.call(paste0, expand.grid(vars, ranges)), ".zip")
  zip <- file.path(outdir, basename(url))

  for (i in 1:12){
    download.file(url[i], zip[i])
  }

  # get all the zip files
  zipF <- list.files(path = outdir, pattern = "*.zip", full.names = TRUE)
  # unzip all your files
  ldply(.data = zipF, .fun = unzip, exdir = outdir)
  #remove all datapoints before landsat (1985)
  file.remove(list.files(path = outdir, pattern = "1980|1981|1982|1983|1984|2019", full.names = TRUE))
}

mk_list <- function(x)  {as.list(intersect(list.files(path = outdir, pattern = "*.tif", full.names = TRUE),
                                           list.files(path = outdir, pattern = x, full.names = TRUE)))}

vars <- c("tmax", "tmin")
for (var in vars){
  list_clim <- mk_list(var)
  # read the files
  raster_clim <- lapply(list_clim, rast)

  yearly <- c()
  for (i in seq(1, 408, 12)){
    tst <- raster_clim[[i]]
    print(i)
    for (j in (i+1):(i+11)){
      tst <- tst + raster_clim[[j]]
      print(j)
    }
    yearly <- c(yearly, tst)
  }
  yearly <- rast(yearly)
  cropped <- crop(yearly, regrowth_mask)
  raster_clim <- resample(cropped, regrowth_mask, method='near')
  raster_clim_masked <- mask(raster_clim, regrowth_mask)
  writeRaster(raster_clim_masked, filename=paste0(var, '_', location, '_santoro.tif'))
}

#resolution is about 4.5 km.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  Soil ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SNUM - Soil mapping unit number

soil <- sf::st_read(dsn = paste0(getwd(),"/soil/") , layer="DSMW")

brazil_soil = soil[soil$COUNTRY == "BRAZIL"]

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





