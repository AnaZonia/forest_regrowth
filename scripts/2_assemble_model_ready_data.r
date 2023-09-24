library(terra) # handling spatial data
library(tidyverse)
setwd("/home/aavila/forest_regrowth")

######################################################################
#################      unify data into central df   ##################
######################################################################
fire <- rast('./data/model_ready_rasters/0000000000-0000095232_fire_history_santoro.tif')
last_LU <- rast('./data/last_LU_0000000000-0000095232_santoro.tif')

tst <- rast('./data/santoro_regrowth.tif')

regrowth_paper <- rast('./data/secondary_forest_age_v2_2018/secondary_forest_age_v2_2018-0000000000-0000065536.tif')
regrowth <- crop(regrowth_paper, tst)
regrowth[regrowth == 0] <- NA

santoro_raster <- rast('./data/santoro/N00W050_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2018-fv4.0.tif')
regrowth_paper <- crop(regrowth_paper, santoro_raster)
santoro_raster <- resample(santoro_raster, regrowth,'near')

santoro_sd <- rast('./data/santoro/N00W050_ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-2018-fv4.0.tif')
santoro_sd <- resample(santoro_sd, regrowth,'near')

all_data_santoro <- c(santoro_raster, regrowth, santoro_sd)

tst <- rast_to_df(all_data_santoro)
colnames(tst) <- c('agbd', 'age', 'sd')

#santoro_raster <- mask(santoro_raster, regrowth)
#tst <- readRDS('santoro_ESA_alldata.rds')

# getting percent of mature forest cover within x neighboring patches:
mature_mask <- rast('./mapbiomas/mature_masks/0000000000-0000095232_mature_mask.tif')
mature_mask <- terra::crop(mature_mask, regrowth)
mature_mask[mature_mask > 0] <- 1
mature_mask <- subst(mature_mask, NA, 0)

range <- 21
mature_sum <- focal(mature_mask, range, sum, na.rm = TRUE) # run focal on these areas
mature_sum <- mask(mature_sum, regrowth) # select only the sum of areas surrounding secondary patches.

### save all data as a unified dataframe for easy visualization and modelling

all_data_santoro <- c(santoro_raster, regrowth, total_prec, total_temp, fire, last_LU)

file_name <- paste0(tempfile(), "_.tif")
# the stack is too big to convert right away
lu_tile <- makeTiles(all_data, c(2000,2000), file_name, na.rm = TRUE, overwrite=TRUE) 
lu_tile <- lapply(lu_tile, rast)

rast_to_df <- function(raster){
  coords <- crds(raster, df=FALSE, na.rm=TRUE)
  values_stack <- terra::extract(raster, coords, cells=FALSE, method="simple")
  central_df <- values_stack[complete.cases(values_stack), ]
  return(central_df)
}

data <- lapply(lu_tile, rast_to_df)
data <- bind_rows(data)

#############################

santoro_raster[santoro_raster == 0] <- NA

rasters <- c(regrowth, santoro_raster, aet, cwd)

common_extent <- ext(rasters[[1]])
for (i in 2:length(rasters)) {
  common_extent <- intersect(common_extent, ext(rasters[[i]]))
}

cropped_rasters <- lapply(rasters, crop, common_extent)


all_data_santoro <- rast(cropped_rasters)
names(all_data_santoro) <- c('agbd', 'age','aet', 'cwd')
all_data_santoro <- trim(all_data_santoro)

file_name <- paste0(tempfile(), "_.tif")
# the stack is too big to convert right away
lu_tile <- makeTiles(all_data_santoro, c(2000,2000), file_name, na.rm = TRUE, overwrite=TRUE) 
lu_tile <- lapply(lu_tile, rast)

rast_to_df <- function(raster){
  coords <- crds(raster, df=FALSE, na.rm=TRUE)
  values_stack <- terra::extract(raster, coords, cells=FALSE, method="simple")
  central_df <- values_stack[complete.cases(values_stack), ]
  return(central_df)
}

data <- lapply(lu_tile, rast_to_df)
data_raw <- bind_rows(data)
#colnames(data) <- c('agbd', 'age', 'prec', 'temp', 'total_fires', 'ts_fire', 'last_LU')

# The issue with matching UTM coordinates directly is a lot of points are lost. what we are plotting is actually a very small percentage of what is actually matching.
# with raster::resample, we can use nearest-neighbor.
# Making rasters from MAPBIOMAS data
lulc <- rast('./model_ready_rasters/0000000000-0000095232_lulc_history.tif') 
regrowth <- rast('./model_ready_rasters/0000000000-0000095232_forest_age.tif')
fire <- rast('./model_ready_rasters/0000000000-0000095232_fire_history.tif')

# Making rasters from GEDI and soil data
# GEDI and soil data is irregular and can't be converted directly into a regular raster.
# making a raster from irregular data, using another raster as reference of size and resolution.
# df must be in form [lon;lat;data].

GEDI <- readRDS(paste0('./dataframes/',"0000000000-0000095232_GEDI.rds"))
  GEDI <- GEDI[,c(3, 2, 5)]
  colnames(GEDI) <- c('lon', 'lat', 'agbd')

soil <- readRDS('./soil/soil.rds') #loads in country-wide soil data
  colnames(soil) <- c('lon', 'lat', 'type')
  soil$type <- as.factor(soil$type)
  # since soil is categorical:
  soil$numtype <- as.numeric(soil$type)

proj <- "+proj=longlat +elips=WGS84"
GEDI_vect <- terra::vect(GEDI[,c("lon", "lat", "agbd")], crs = proj)
GEDI_raster <- terra::rasterize(GEDI_vect, regrowth, field = "agbd")

soil_vect <- terra::vect(soil[,c("lon", "lat", "numtype")], crs = proj)
soil_raster <- terra::rasterize(soil_vect, GEDI_raster, field = "numtype")



