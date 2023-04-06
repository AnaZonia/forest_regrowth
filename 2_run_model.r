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

######################################################################
#################      unify data into central df   ##################
######################################################################
# The issue with matching UTM coordinates directly is a lot of points are lost. what we are plotting is actually a very small percentage of what is actually matching.
# with raster::resample, we can use nearest-neighbor.
# Making rasters from MAPBIOMAS data
lulc <- rast('./model_ready_rasters/0000000000-0000095232_lulc_history.tif') 
regrowth <- rast('./model_ready_rasters/0000000000-0000095232_forest_age.tif')
fire <- rast('./model_ready_rasters/0000000000-0000095232_fire_historu.tif')

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
soil_rasample <- resample(soil_raster, regrowth_raster, method = 'near')

temp <- rast(paste0('./worldclim_dataframes/','temp_BRA_mean.tif'))
prec <- rast(paste0('./worldclim_dataframes/','prec_BRA_mean.tif'))
prec_sd <- rast(paste0('./worldclim_dataframes/','prec_BRA_sd.tif'))
temp_sd <- rast(paste0('./worldclim_dataframes/','temp_BRA_sd.tif'))

# resampling the rasters
# soil and GEDI don't need to be resampled as they have already been rasterized
# with regrowth_raster as reference
prec_resampled <- resample(prec,regrowth_raster,method = 'near')
temp_resampled <- resample(temp,regrowth_raster,method = 'near')
prec_sd_resampled <- resample(prec_sd,regrowth_raster,method = 'near')
temp_sd_resampled <- resample(temp_sd,regrowth_raster,method = 'near')
fire_resampled <- resample(fire_raster,regrowth_raster,method = 'near')

GEDI_resampled <- resample(GEDI_raster,regrowth_raster,method = 'near')
lulc_resampled <- resample(lulc_raster,regrowth_raster,method = 'near')

lulc_masked <- mask(lulc_raster, regrowth_raster, maskvalue = NaN, updatevalue = NA)
GEDI_masked <- mask(GEDI_resampled, regrowth_raster, maskvalue = NaN, updatevalue = NA)
fire_masked <- mask(fire_resampled, regrowth_raster, maskvalue = NaN, updatevalue = NA)

global(regrowth_raster, 'notNA')
global(lulc_masked$soy, 'notNA')
global(GEDI_masked, 'notNA')


sum_prec <- mask(sum(prec_resampled), regrowth_raster)
sum_temp <- mask(sum(temp_resampled), regrowth_raster)
sum_prec_sd <- mask(sum(prec_sd_resampled), regrowth_raster)
sum_temp_sd <- mask(sum(temp_sd_resampled), regrowth_raster)

names(GEDI_raster) <- 'agbd'
names(regrowth_raster) <- 'age'
names(fire_resampled) <- c('num_fires', 'ts_last_fire')
names(soil_raster) <- 'soil_type'
names(sum_prec) <- 'sum_prec'
names(sum_temp) <- 'sum_temp'
names(sum_prec_sd) <- 'sum_prec_sd'
names(sum_temp_sd) <- 'sum_temp_sd'

######################################################################
#################        Sampling     ##################
######################################################################

agbd <- terra::as.data.frame(GEDI_resampled, xy=TRUE, na.rm=TRUE)
regrowth <- terra::as.data.frame(regrowth_raster, xy=TRUE, cells = TRUE)
#lulc <- terra::as.data.frame(lulc_resampled, xy=TRUE, cells = TRUE)
lulc_masked_df <- terra::as.data.frame(lulc_masked, xy=TRUE, cells = TRUE)
fire <- terra::as.data.frame(fire_masked, xy=TRUE, cells = TRUE)




mature_mask <- rast('./mapbiomas/mature_masks/0000000000-0000095232_mask.tif')
mature_mask <- terra::crop(mature_mask, regrowth_raster)
# nrow <- 10185
# ncol <- 16055
mature_mask[mature_mask > 0] <- 1
mature_mask <- subst(mature_mask, NA, 0)

regrowth_subset <- lapply(split(regrowth, regrowth$age),
   function(subdf) subdf[sample(1:nrow(subdf), 3),]
)

regrowth_subset <- do.call('rbind', regrowth_subset)
head(regrowth_subset)

# look for indices in the mature_mask
# 
# convert straight directions to x and y values
covered_area <- function(reg_index, b){
  #b <- 40
  #reg_index <- regrowth_subset$cell[1]
  col_num <- (reg_index - (reg_index %% nrow(mature_mask)))/nrow(mature_mask) + 1
  row_num <- reg_index %% nrow(mature_mask)
  x <- mature_mask[(row_num-b):(row_num+b), (col_num-b):(col_num+b), drop=FALSE]
  return(global(x, sum))
}

sums <- lapply(regrowth_subset$cell, covered_area, 400)
# converting the cells from one to the other is not needed here since I am using two rasters
# add total of b neighboring cells

# 41% of cells in mature_mask are not NA
global(lulc_raster$pasture, 'isNA')
global(regrowth_raster, 'isNA')


######################################################################
#################        passing into the model     ##################
######################################################################


# fit the sum of total temperature/tmp/Rtmp6q4gQ6/vscode-R/plot.png
# fit soil as categorical


# Information within lulc_raster:
# lulc_raster[[1]] <- pasture
# this file has no values for sugar and coffee plantations, so this area is only considering 4 land use types.
# lulc_raster[[1]] <- pasture
# lulc_raster[[2]] <- soy
# lulc_raster[[5]] <- other_perennial
# lulc_raster[[6]] <- other_annual
# lulc_raster[[7]] <- time since pasture
# lulc_raster[[8]] <- time since soy
# lulc_raster[[9]] <- time since other_perennial
# lulc_raster[[10]] <- time since other_annual

# test_value[[1]] <- agbd
# test_value[[2]] <- age
# test_value[[3]] <- num_fires
# test_value[[4]] <- ts_last_fire
# test_value[[5]] <- soil_type
# test_value[[6]] <- sum_prec
# test_value[[7]] <- sum_temp
# test_value[[8]] <- sum_prec_sd
# test_value[[9]] <- sum_temp_sd
# test_value[[10]] <- pasture
# test_value[[11]] <- soy
# test_value[[12]] <- coffee
# test_value[[13]] <- sugar
# test_value[[14]] <- other_perennial
# test_value[[15]] <- other_annual
# test_value[[16]] <- time since pasture
# test_value[[17]] <- time since soy
# test_value[[18]] <- time since other_perennial
# test_value[[19]] <- time since other_annual


G <- function(pars) {
  # Extract parameters of the model
  Gmax = pars[1] #asymptote
  k = pars[10]*test_value[[3]] ^ (pars[11]*test_value[[4]]) # num_fires ^ ts_last_fire
  +pars[2]*test_value[[10]] ^ (pars[3] * test_value[[16]]) # pasture ^ ts_pasture
  +pars[4]*test_value[[11]] ^ (pars[5] * test_value[[17]]) # soy ^ ts_soy
  +pars[6]*test_value[[14]] ^ (pars[7] * test_value[[18]]) # other_perennial ^ ts_other_perennial
  +pars[8]*test_value[[15]] ^ (pars[9] * test_value[[19]]) # other_annual ^ ts_other_annual
  + exp(-pars[12]*test_value[[6]]) + pars[13]*test_value[[7]] # precipitation ; temperature
  # Prediction of the model
  Gmax * (1 - exp(-k*(test_value[[2]])))
}

####################
#finding the right parameters
par0 = c(50,  0.1, 0.1, 0.1,  0.1,
         0.1, 0.1, 0.1,  0.1,  0.1,
         0.1, 0.1, 0.1,  0.1, 0.1)
Gpred <- G(par0)
Gpred

NLL = function(pars) {
  # Values prediced by the model
  pars <- par0
  if(pars[length(pars)] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  Gpred = G(pars)
  #print(Gpred)
  # Negative log-likelihood 
  fun = -sum(dnorm(x = as.vector(test_value[[1]]), mean = as.vector(Gpred), sd = pars[length(pars)], log = TRUE), na.rm = TRUE)
  #print(pars)
  return(fun)
}

NLL(par0)

o = optim(par = par0, fn = NLL, control = list(parscale = abs(par0)), 
           hessian = FALSE, method = "BFGS")
print(o)

meth0 = c("Nelder-Mead", "BFGS", "CG")
for (i in meth0){
  o = optim(par = par0, fn = NLL, control = list(parscale = abs(par0)), 
             hessian = FALSE, method = i)
  print(i)
  print(o)
}

pred = G(o$par[1:13])

plot(central_df$agbd, pred, abline(0,1))

