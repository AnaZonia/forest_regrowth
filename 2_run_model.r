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
  # since soil is categorical:
  soil$numtype <- as.numeric(soil$type)

prec <- rast('./model_ready_rasters/prec_0000000000-0000095232.tif')
#tmin <- rast('./model_ready_rasters/tmin_0000000000-0000095232.tif')
#tmax <- rast('./model_ready_rasters/tmax_0000000000-0000095232.tif')
#temp <- tmax-tmin
#writeRaster(temp, paste0('./model_ready_rasters/temp_0000000000-0000095232.tif'))
temp <- rast('./model_ready_rasters/temp_0000000000-0000095232.tif')

proj <- "+proj=longlat +elips=WGS84"
GEDI_vect <- terra::vect(GEDI[,c("lon", "lat", "agbd")], crs = proj)
GEDI_raster <- terra::rasterize(GEDI_vect, regrowth, field = "agbd")
soil_vect <- terra::vect(soil[,c("lon", "lat", "numtype")], crs = proj)
soil_raster <- terra::rasterize(soil_vect, regrowth, field = "numtype")

fire_cropped <- crop(fire, GEDI_raster)
GEDI_raster_cropped <- crop(GEDI_raster, fire)
regrowth_cropped <- crop(regrowth, GEDI_raster_cropped)
lulc_cropped <- crop(lulc, GEDI_raster_cropped)
soil_cropped <- crop(soil_raster, GEDI_raster_cropped)
temp_cropped <- crop(temp, GEDI_raster_cropped)
prec_cropped <- crop(prec, GEDI_raster_cropped)

# names(GEDI_raster) <- 'agbd'
# names(regrowth_cropped) <- 'age'
# names(fire_resampled) <- c('num_fires', 'ts_last_fire')
# names(soil_raster) <- 'soil_type'
# names(sum_prec) <- 'sum_prec'
# names(sum_temp) <- 'sum_temp'
# names(sum_prec_sd) <- 'sum_prec_sd'
# names(sum_temp_sd) <- 'sum_temp_sd'

# to save time as I'm figuring things out

all_data_raster <- c(GEDI_raster_cropped, regrowth_cropped, fire_cropped, lulc_cropped, soil_cropped, temp_cropped, prec_cropped)

writeRaster(all_data_raster, 'all_data_raster.tif')

all_data_raster <- rast('all_data_raster.tif')
regrowth_cropped <- all_data_raster[[2]]
regrowth_cropped
######################################################################
#################        Sampling     ##################
######################################################################
total_temp <- app(temp_cropped, sum)
total_prec <- app(prec_cropped, sum)
total_years <- lulc_cropped[[6:9]]
last_LU <- lulc_cropped[[1]]

# getting percent of mature forest cover within x neighboring patches:
# would like to get values within 300m radius (range of dispersal - cite)
mature_mask <- rast('./mapbiomas/mature_masks/0000000000-0000095232_mature_mask.tif')
mature_mask <- terra::crop(mature_mask, regrowth_cropped)
mature_mask[mature_mask > 0] <- 1
mature_mask <- subst(mature_mask, NA, 0)

window_size <- 21 # (must be an odd number)

mature_sum <- focal(mature_mask, window_size, sum, na.rm = TRUE)
mature_sum <- mask(mature_sum, regrowth_cropped)
# think again of how to get the asymptote
# think of how to get the percent cover of nearby mature forest
# create a mask with the region around the secondary forest patches.

# all other values will be deleted. focal should run faster
# mask with only the cells around the vicinity of a secondary forest

# create vectors with the coordinates around each patch
# extract said values
# add values to the corresponding coordinate in a blank new raster.
secondary_patches <- cells(regrowth_cropped)

secondary_patches <- 13 #cells(regrowth_cropped)
remainder_j <- secondary_patches %% ncol(regrowth_cropped)
multiple_i <- (secondary_patches - remainder_j)/ncol(regrowth_cropped)
tst <- data.frame(secondary_patches, remainder_j, multiple_i)
mat_total <- 0
for (k in -range:range){
  min_indx <- (multiple_i+k)*ncol(mature_mask)+(remainder_j-range)
  max_indx <- ((multiple_i+k)*ncol(mature_mask)+(remainder_j+range))
  mat_total <- rm_total + sum( mature_mask[min_indx : max_indx])
}

tst <- vect(secondary_patches)

rk <- matrix(1, 101, 101)
tst <- adjacent(regrowth_cropped, cells=secondary_patches, directions=rk) # why does bad alloc

# this should be faster than the focal function as it only operates on areas that have secondary forest.
function <- total_mature_patches(regrowth_cropped, mature_mask, range){  # range is the distance from the secondary forest
  range <- 21
  m <- matrix(1:25,5,5)
  rm <- rast(m)


}



  m <- matrix(1:36,6,6)
  rm <- rast(m)

  m2 <- matrix(0,6,6)
  rm2 <- rast(m2)
  rm2[15] <- 1
  rm2[22] <- 1

  secondary_patches <- c(15, 22) #cells(regrowth_cropped)
  remainder_j <- secondary_patches %% ncol(rm)
  multiple_i <- (secondary_patches - remainder_j)/ncol(rm)
  tst <- data.frame(secondary_patches, remainder_j, multiple_i)

  range <- 2
  
  apply(tst, 1, function(multiple_i, remainder_j) {
     mat_total <- 0
    for (k in -range:range){
      min_indx <- (multiple_i+k)*ncol(mature_mask)+(remainder_j-range)
      max_indx <- ((multiple_i+k)*ncol(mature_mask)+(remainder_j+range))
      mat_total <- rm_total + sum( mature_mask[min_indx : max_indx])
    }    return(mat_total)
    }
  )



# mean biomass of surrounding mature forests around 3km (~100 pixels of 30m)
surrounding_forest <- 101 # window size
mature_biomass <- mask(GEDI_raster_cropped, mature_mask)
mature_mean_biomass <- focal(mature_biomass, surrounding_forest, mean, na.rm = TRUE)
# biomass_surr_mat <- mask(mature_mean_biomass, regrowth_cropped) # mean biomass of surrounding mature forests for every secondary forest patch

all_data <- c(GEDI_raster_cropped, regrowth_cropped, fire_cropped, total_years, last_LU, total_temp, total_prec)

coords <- crds(all_data[[1]], df=FALSE, na.rm=TRUE)
values_stack <- terra::extract(all_data, coords, cells=FALSE, method="simple")
central_df <- values_stack[complete.cases(values_stack), ]
colnames(central_df)[1:4] <- c('agbd', 'age', 'total_fires', 'ts_fire')
colnames(central_df)[10:11] <- c('temp', 'prec')

central_df$Bmax <- 400
#saveRDS(central_df, 'central_df.rds')

# normalize central_df
minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

normalized_df <- as.data.frame(lapply(central_df, minMax))
normalized_df$Bmax <- 400

######################################################################
#################        passing into the model     ##################
######################################################################

# create another list that sets to zero if things are not present

G <- function(pars) {
  # Extract parameters of the model
  E = pars[1]*normalized_df$prec + pars[2]*normalized_df$temp #+ pars[3]*central_df$mature_sum
  LU = normalized_df$total_fires * exp(pars[4]*normalized_df$ts_fire + pars[5]*normalized_df$last_LU)
  k = E + LU
  # Prediction of the model
  return ( normalized_df$Bmax * (1 - exp(-k)) )
}

par0 = c(0.5,0.5,0.5,0.5,0.5)
Gpred <- G(par0)
Gpred

LU = normalized_df$total_fires * exp(0.05*normalized_df$ts_fire + 0.05*normalized_df$last_LU)



# k is going to be very large so the term is going very tiny.
# walk through the function. break it down into component parts.
# anything below 20 is going to be tiny.

# asymptote is being repeated.
# have so that the different selections have variance.
# normalize the variables so they dont vary wildly. divide the parameter values.

####################
#finding the right parameters

NLL = function(pars) {
  # Values prediced by the model
  pars <- par0
  if(pars[length(pars)] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  # before returning likelihood value, check if it is NA or inf or something problematic
  Gpred = G(pars)
  #print(Gpred)
  # Negative log-likelihood 
  fun = -sum(dnorm(x = central_df$agbd - Gpred, mean = 0, sd = pars[length(pars)], log = TRUE), na.rm = TRUE)
  return(fun)
}

# are the likelihoods varying?
# if the surface 
NLL(par0)

o = optim(par = par0, fn = NLL, hessian = FALSE, method = "BFGS")
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

######################################





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
