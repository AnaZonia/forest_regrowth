####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the district of Paragominas, in the Western Amazon.
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################

library(terra) # handling spatial data
library(data.table) #to use data.table objects in the rasterFromXYZ_irr() function
library(pbapply) #progress bar for apply family of functions
library(sf)
setwd("/home/aavila/forest_regrowth")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")

##################      Switches     ##################

import_mapbiomas = TRUE

######################################################################
#################      unify data into central df   ##################
######################################################################
# The issue with matching UTM coordinates directly is a lot of points are lost. what we are plotting is actually a very small percentage of what is actually matching.
# with raster::resample, we can use nearest-neighbor.

# Making rasters from MAPBIOMAS data

lulc <- readRDS(paste0('./dataframes/', '0000000000-0000095232_lulc_history.rds'))
regrowth <- readRDS(paste0('./dataframes/', '0000000000-0000095232_regrowth_history.rds'))
fire <- readRDS(paste0('./dataframes/','0000000000-0000095232_fire_history.rds'))

lulc_raster <- rast(lulc[,c(1:12)], type="xyz") #lon, lat, values    , crs = "+init=epsg:4326"
regrowth_raster <- rast(regrowth[,c(1,2,4)], type="xyz")
fire_raster <- rast(fire[,c(1:4)], type="xyz")



# Making rasters from GEDI and soil data
# GEDI and soil data is irregular and can't be converted directly into a regular raster.
# making a raster from irregular data, using another raster as reference of size and resolution.
# df must be in form [lon;lat;data].

GEDI <- readRDS(paste0('./dataframes/',"0000000000-0000095232_GEDI.rds")) #for some reason, the method doesn't work with this file?
  GEDI <- GEDI[,c(3, 2, 5)]
  colnames(GEDI) <- c('lon', 'lat', 'agbd')

soil <- readRDS('./soil/soil.rds') #loads in country-wide soil data
  colnames(soil) <- c('lon', 'lat', 'type')
  soil$type <- as.factor(soil$type)

# I am having a hard time making this a single function - to be discovered why.
proj <- "+proj=longlat +elips=WGS84"
GEDI_vect <- terra::vect(GEDI[,c("lon", "lat", "agbd")], crs = proj)
GEDI_raster <- terra::rasterize(GEDI_vect, regrowth_raster, field = "agbd")

soil_vect <- terra::vect(soil[,c("lon", "lat", "type")], crs = proj)
soil_raster <- terra::rasterize(soil_vect, regrowth_raster, field = "type")

temp <- rast(paste0('./worldclim_dataframes/','temp_BRA_mean.tif'))
prec <- rast(paste0('./worldclim_dataframes/','prec_BRA_mean.tif'))

# resampling the rasters
# soil and GEDI don't need to be resampled as they have already been rasterized
# with regrowth_raster as reference
prec_resampled <- resample(prec,regrowth_raster,method = 'near')
temp_resampled <- resample(prec,regrowth_raster,method = 'near')
fire_resampled <- resample(fire_raster,regrowth_raster,method = 'near')

names(GEDI_raster) <- 'agbd'
names(regrowth_raster) <- 'age'
names(fire_resampled) <- c('num_fires', 'ts_last_fire')
names(prec_resampled)
names(temp_resampled)
names(soil_raster)

central_df <- c(GEDI_raster, regrowth_raster, lulc_resampled,
                fire_resampled, prec_resampled, temp_resampled, soil_raster)

######################################################################
#################        passing into the model     ##################
######################################################################

# fit the sum of total temperature
# fit soil as categorical


G <- function(pars, columns) {
  # Extract parameters of the model
  Gmax = pars[1] #asymptote
  for (i in 1:length(columns)){
    k = pars[2]*other_perennial+pars[3]*pasture+pars[4]*soy+pars[5]*other_annual + pars[6]*other_annual + pars[7]*num_fires^(pars[8]*last_burn)
  }
  # Prediction of the model
  Gmax * (1 - exp(-k*(forest_age)))
}

####################
#finding the right parameters
par0 = c(50, 0.0001,  0.0001,  0.001, 0.001,  0.001, 0.1, 0.1)

Gpred = G(par0, central_df$other_perennial, central_df$other_annual, central_df$pasture, central_df$soy, central_df$forest_age, central_df$num_fires, central_df$last_burn)

# dat$tavg, dat$prec, 

NLL = function(pars, dat) {
  # Values prediced by the model
  if(pars[9] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  Gpred = G(pars, dat$other_perennial, dat$other_annual, dat$pasture, dat$soy, dat$forest_age, dat$num_fires, dat$last_burn)
  #print(Gpred)
  # Negative log-likelihood 
  fun = -sum(dnorm(x = dat$agbd, mean = Gpred, sd = pars[9], log = TRUE))
  #print(pars)
  return(fun)
}

par0 = c(50, 0.001, 0.0001,  0.0001,  0.0001,  0.0001,  0.0001, 0.1, 0.1)
NLL(par0, dat=central_df)

o = optim(par = par0, fn = NLL, dat = central_df, control = list(parscale = abs(par0)), 
           hessian = FALSE, method = "CG")
print(o)

meth0 = c("Nelder-Mead", "BFGS", "CG")
for (i in meth0){
  o = optim(par = par0, fn = NLL, dat = central_df, control = list(parscale = abs(par0)), 
             hessian = FALSE, method = i)
  print(i)
  print(o)
}

#  central_df$tavg, central_df$prec, 

pred = G(o$par[1:9], central_df$other_perennial, central_df$other_annual, central_df$pasture, central_df$soy, central_df$forest_age, central_df$num_fires, central_df$last_burn)

plot(central_df$agbd, pred, abline(0,1))

