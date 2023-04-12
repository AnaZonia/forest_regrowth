####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the district of Paragominas, in the Western Amazon.
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################

library(terra) # handling spatial data
library(data.table) #to use data.table objects in the rasterFromXYZ_irr() function
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

temp <- rast(paste0('./model_ready_rasters/','temp_BRA_mean.tif'))
prec <- rast(paste0('./model_ready_rasters/','prec_BRA_mean.tif'))
prec_sd <- rast(paste0('./model_ready_rasters/','prec_BRA_sd.tif'))
temp_sd <- rast(paste0('./model_ready_rasters/','temp_BRA_sd.tif'))

# just for testing sake, for now:
tst_raster <- rast('./dataframes/0000000000-0000095232_regrowth_raster.tif')
proj <- "+proj=longlat +elips=WGS84"
GEDI_vect <- terra::vect(GEDI[,c("lon", "lat", "agbd")], crs = proj)
GEDI_raster <- terra::rasterize(GEDI_vect, tst_raster, field = "agbd")
soil_vect <- terra::vect(soil[,c("lon", "lat", "type")], crs = proj)
soil_raster <- terra::rasterize(soil_vect, tst_raster, field = "type")

# x = 1.5, y = 1.9
# x = 1.6 ...
# as.integer, divide by cell size of landsat, 1.6 goes to 1.5
# gives cell index. relate everything to cell indices.
# 

fire_cropped <- crop(fire, GEDI_raster)
GEDI_raster_cropped <- crop(GEDI_raster, fire)
regrowth_cropped <- crop(regrowth, GEDI_raster_cropped)
lulc_cropped <- crop(lulc, GEDI_raster_cropped)
soil_cropped <- crop(soil_raster, GEDI_raster_cropped)

# resampling the rasters
prec_resampled <- resample(prec,regrowth_cropped,method = 'near')
temp_resampled <- resample(temp,regrowth_cropped,method = 'near')
prec_sd_resampled <- resample(prec_sd,regrowth_cropped,method = 'near')
temp_sd_resampled <- resample(temp_sd,regrowth_cropped,method = 'near')
# isn't working perfectly. some areas with land use don't have climate data.

sum_prec <- mask(sum(prec_resampled), regrowth_cropped)
sum_temp <- mask(sum(temp_resampled), regrowth_cropped)
sum_prec_sd <- mask(sum(prec_sd_resampled), regrowth_cropped)
sum_temp_sd <- mask(sum(temp_sd_resampled), regrowth_cropped)

# names(GEDI_raster) <- 'agbd'
# names(regrowth_cropped) <- 'age'
# names(fire_resampled) <- c('num_fires', 'ts_last_fire')
# names(soil_raster) <- 'soil_type'
# names(sum_prec) <- 'sum_prec'
# names(sum_temp) <- 'sum_temp'
# names(sum_prec_sd) <- 'sum_prec_sd'
# names(sum_temp_sd) <- 'sum_temp_sd'

######################################################################
#################        Sampling     ##################
######################################################################

total_years <- lulc_cropped[[6:9]]
last_LU <- lulc_cropped[[1]]
total_fires <- fire_cropped[[1]]

tst <- c(GEDI_raster_cropped, regrowth_cropped, total_fires, total_years, last_LU)
tst <- mask(tst, GEDI_raster_cropped)
tst <- mask(tst, regrowth_cropped)

coords <- crds(tst[[1]], df=FALSE, na.rm=TRUE)
values_stack <- terra::extract(tst, coords, cells=FALSE, method="simple")
central_df <- values_stack[complete.cases(values_stack), ]
colnames(central_df)[1:3] <- c('agbd', 'age', 'total_fires')

mature_mask <- rast('./mapbiomas/mature_masks/0000000000-0000095232_mature_mask.tif')
mature_mask <- terra::crop(mature_mask, regrowth_cropped)
mature_mask[mature_mask > 0] <- 1
mature_mask <- subst(mature_mask, NA, 0)

regrowth_subset <- lapply(split(central_df, central_df$age),
   function(subdf) subdf[sample(1:nrow(subdf), 3),]
)
# by (central_df, central_df$age)
# tapply 1,nrow(centraldf), central_df$age
# get a vector of rownumbers
# regrowth_subset, if pos = tapply 1,nrow(centraldf), central_df$age, 
# regrowth_subset [pos]

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

######################################################################
#################        passing into the model     ##################
######################################################################


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

