####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the district of Paragominas, in the Western Amazon.
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################

library(terra) # handling spatial data
library(raster)
library(data.table) #to use data.table objects in the rasterFromXYZ_irr() function
library(pbapply) #progress bar for apply family of functions

setwd("/home/aavila/forest_regrowth/dataframes")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")

######################################################################
#################      unify data into central df   ##################
######################################################################

# The issue with matching UTM coordinates directly is a lot of points are lost. what we are plotting is actually a very small percentage of what is actually matching.
# with raster::resample, we can use nearest-neighbor.

lulc <- readRDS('0000000000-0000095232_lulc_history.rds')
age <- readRDS('0000000000-0000095232_forest_age.rds')
burn <- readRDS('0000000000-0000095232_burn_history.rds')

# GEDI_mid_amazon <- readRDS('GEDI_midamazon_dfunified.rds')
#   GEDI_mid_amazon <- GEDI_mid_amazon[,c(3,2,1)]
GEDI <- readRDS("0000000000-0000095232_GEDI.rds") #for some reason, the method doesn't work with this file?
  GEDI <- GEDI[,c(3,2,5,6,7,8)]
  colnames(GEDI) <- c('lon', 'lat', 'agbd', 'zone', 'x', 'y')
  GEDI <- setDT(GEDI)

lons <- c(GEDI[,1], age[,1], burn[,1], lulc[,1])
lats <- c(GEDI[,2], age[,2], burn[,2], lulc[,2])
e_min <- c(max(lapply(lons, min)), min(lapply(lons, max)), max(lapply(lats, min)), min(lapply(lats, max)))
e_min # extent object containing the smallest possible extent encompassing all dataframes.
# crop the dataframes with extent different from e_min

dfs <- c(GEDI, age, burn, lulc)

tst = lapply(lons, min)
min(tst)

range(lulc$lon)
range(lulc$lat)
range(age$lon)
range(age$lat)
range(burn$lon)
range(burn$lat)
range(GEDI$lon)
range(GEDI$lat)

# temp_raster <- rasterFromXYZ(temp[,c(1:37)], crs = "+init=epsg:4326")
# prec_raster <- rasterFromXYZ(prec[,c(1:37)], crs = "+init=epsg:4326")
lulc_raster <- rasterFromXYZ(lulc[,c(1:12)], crs = "+init=epsg:4326")
age_raster <- rasterFromXYZ(age[,c(1,2,4)], crs = "+init=epsg:4326")
burn_raster <- rasterFromXYZ(burn[,c(1:4)], crs = "+init=epsg:4326")

# checking extents for equivalency - finding the extent that will encompass all.
# due to variation within the data, there's some risk of small differences.
# extent(temp_raster)
# extent(lulc_raster)
# extent(age_raster)
# extent(burn_raster)
raster_list <- c(lulc_raster, burn_raster, age_raster, temp_resampled, prec_resampled)
raster_list <- pbapply::pblapply(raster_list, terra::crop, e_min)
writeRaster(lulc_raster, '0000000000-0000095232_lulc_raster.tif')
writeRaster(burn_raster, '0000000000-0000095232_burn_raster.tif')
writeRaster(age_raster, '0000000000-0000095232_age_raster.tif')
writeRaster(prec_resampled, '0000000000-0000095232_prec_resampled.tif')
writeRaster(temp_resampled, '0000000000-0000095232_temp_resampled.tif')


soil <- readRDS('./soil/soil.rds')
  colnames(soil) <- c('lon', 'lat', 'type')
  soil <- soil[max(age$lon) > soil$lon, ] # this code needs to be simplified
  soil <- soil[min(age$lon) < soil$lon, ]
  soil <- soil[max(age$lat) > soil$lat, ]
  soil <- soil[min(age$lat) < soil$lat, ]
  soil <- setDT(soil)

# preds <- c(temp[,c(1:37)], prec[,c(1:37)], age[,c(1,2,4)])
# preds2 <- pbapply::pblapply(preds, rasterFromXYZ, crs = "+init=epsg:4326")  <------ not sure why this is returning an incorrect number of dimensions error?

# resampling the rasters
temp_resampled <- resample(temp_raster,age_raster,method = 'ngb') #nearest neighbor
prec_resampled <- resample(prec_raster,age_raster,method = 'ngb')


# GEDI and soil data is irregular and can't be converted directly into a regular raster.
# making a raster from irregular data, using another raster as reference of size and resolution.
# df must be in form [lon;lat;data].
rasterFromXYZ_irr <- function(df, ref_raster, col){   # dataframe to be converted, reference raster
  df <- GEDI
  ref_raster <- age_raster
  e <- extent(min(df$lon), max(df$lon), min(df$lat), max(df$lat))
  # set up an 'empty' raster
  r <- raster(e, ncol = ncol(ref_raster), nrow = nrow(ref_raster))
  df[,3] <- as.factor(df[,3])
  x <- rasterize(df[, 1:2], r, col)
  return(x)
}

GEDI_raster <- rasterFromXYZ_irr(GEDI, age_raster, GEDI$agbd)
soil_raster <- rasterFromXYZ_irr(soil, age_raster, soil$type)
#writeRaster(GEDI_raster, '0000000000-0000095232_GEDI_tst.tif')
GEDI_raster <- raster('0000000000-0000095232_GEDI_tst.tif')

writeRaster(central_stack, '0000000000-0000095232_central_stack.tif')

central_stack <- stack(raster_list)
# The other strategy is using the match() function.
# However, this is somehow missing a lot of matches. To be investigated why that is the case.

######################################################################
#################        passing into the model     ##################
######################################################################

# fit the sum of total temperature
# fit soil as categorical
temp_hist <- central_stack
prec_hist <- 
columns <- c(other_perennial, other_annual, pasture, soy,
            ts_other_perennial, ts_other_annual, ts_pasture, ts_soy,
            forest_age, num_fires, last_burn, soil)

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

