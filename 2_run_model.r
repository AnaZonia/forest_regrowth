####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the district of Paragominas, in the Western Amazon.
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################

library(ggplot2)
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
soil_raster <- terra::rasterize(soil_vect, GEDI_raster, field = "numtype")
#soil_raster <- focal(soil_raster, 301, "modal", NAonly=TRUE, na.rm = TRUE)

fire_cropped <- crop(fire, GEDI_raster)
GEDI_raster_cropped <- crop(GEDI_raster, fire)
regrowth_cropped <- crop(regrowth, GEDI_raster_cropped)
lulc_cropped <- crop(lulc, GEDI_raster_cropped)
soil_cropped <- crop(soil_raster, GEDI_raster_cropped) # work on this
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

GEDI_raster <- mask(GEDI_raster, regrowth)
regrowth <- mask(regrowth, GEDI_raster)

all_data <- c(GEDI_raster, regrowth)

coords <- crds(all_data[[1]], df=FALSE, na.rm=TRUE)
values_stack <- terra::extract(all_data, coords, cells=FALSE, method="simple")
agb_forest_age <- values_stack[complete.cases(values_stack), ]

sds <- aggregate(agbd ~ forest_age, agb_forest_age, sd)
means <- aggregate(agbd ~ forest_age, agb_forest_age, mean)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'mean', 'sd')

ggplot(sum_stats,                               # ggplot2 plot with means & standard deviation
       aes(x = age,
           y = median)) + 
  geom_errorbar(aes(ymin = median - sd,
                    ymax = median + sd)) +
  geom_point() + theme(text = element_text(size = 20))  


tst = subset(agb_forest_age, forest_age == 28)
tst = subset(tst, agbd < 25)


######################################################################
#################        Sampling     ##################
######################################################################
total_temp <- app(temp, sum)
total_prec <- app(prec, sum)
total_years <- lulc_cropped[[6:9]]
last_LU <- lulc_cropped[[1]]

# getting percent of mature forest cover within x neighboring patches:
# would like to get values within 300m radius (range of dispersal - cite)
mature_mask <- rast('./mapbiomas/mature_masks/0000000000-0000095232_mature_mask.tif')
mature_mask <- terra::crop(mature_mask, regrowth)
mature_mask[mature_mask > 0] <- 1
mature_mask <- subst(mature_mask, NA, 0)

range <- 21
mature_sum <- focal(mature_mask, range, sum, na.rm = TRUE) # run focal on these areas
mature_sum <- mask(mature_sum, regrowth_cropped) # select only the sum of areas surrounding secondary patches.

# mean biomass of surrounding mature forests around 3km (~100 pixels of 30m)
biomass_range <- 101 # window size
mature_biomass <- mask(GEDI_raster_cropped, mature_mask)  # get only biomass of mature forests
mature_total_biomass <- focal(mature_biomass, biomass_range, sum, na.rm = TRUE)

all_data <- c(GEDI_raster_cropped, regrowth_cropped, fire_cropped, total_years, last_LU, total_temp, total_prec, mature_sum) #, soil_cropped)

coords <- crds(all_data[[1]], df=FALSE, na.rm=TRUE)
values_stack <- terra::extract(all_data, coords, cells=FALSE, method="simple")
central_df <- values_stack[complete.cases(values_stack), ]
colnames(central_df)[1:4] <- c('agbd', 'age', 'total_fires', 'ts_fire')
colnames(central_df)[10:12] <- c('temp', 'prec', 'mat_sum')

#saveRDS(central_df, 'central_df.rds')
central_df <- readRDS('central_df.rds')
central_df <- central_df[central_df$last_LU %in% c(15, 41, 48),]


central_df$last_LU <- factor(central_df$last_LU)
dummy_LU <- as.data.frame(model.matrix(~ central_df$last_LU - 1))
names(dummy_LU) <- c('pasture', 'other_annual', 'other_perennial')

# min-max normalize central_df
minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

normalized_df <- as.data.frame(lapply(central_df[,-c(6,9)], minMax))

#normalized_df <- as.data.frame(lapply(central_df[,-c(6,9)], scale))

normalized_df$Bmax <- 400
head(normalized_df)

data <- cbind(normalized_df, dummy_LU)
head(data)

# plot biomass and rainfall for mature forests

library(mgcv)
fit <- lm(data$agbd ~ data$age)
by(central_df$agbd, central_df$age, summary)
summary(gam(central_df$agbd ~ central_df$age))

jpeg('plot3.jpeg')
plot(data$agbd ~ data$age)
abline(fit)
dev.off()

pred <- predict(fit)
jpeg('pred_obs.jpeg')
plot(pred, data$agbd)
abline(0,1)
dev.off()

##############
regrowth_paper <- rast('./secondary_forest_age_v2_2018/secondary_forest_age_v2_2018-0000000000-0000065536.tif')
regrowth2 <- crop(regrowth_paper, regrowth)
regrowth2[regrowth2 == 0] <- NA

santoro_raster <- rast('N00W050_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2018-fv4.0.tif')
santoro_raster <- crop(santoro_raster, regrowth2)
santoro_raster <- resample(santoro_raster, regrowth2)
santoro_raster <- mask(santoro_raster, regrowth2)
all_data <- c(santoro_raster, regrowth2)

writeRaster(all_data, 'santoro_regrowth.tif')
coords <- crds(all_data[[1]], df=FALSE, na.rm=TRUE)
values_stack <- terra::extract(all_data, coords, cells=FALSE, method="simple")
central_df <- values_stack[complete.cases(values_stack), ]
colnames(central_df)[1:2] <- c('agbd', 'age')

fit <- lm(central_df$agbd ~ central_df$age)
summary(fit)
plot(central_df$agbd ~ central_df$age)
abline(fit)

sds <- aggregate(agbd ~ age, central_df, sd)
means <- aggregate(agbd ~ age, central_df, median)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'median', 'sd')

ggplot(sum_stats,                               # ggplot2 plot with means & standard deviation
       aes(x = age,
           y = median)) + 
  geom_errorbar(aes(ymin = median - sd,
                    ymax = median + sd)) +
  geom_point() + theme(text = element_text(size = 20))  


##################

mature_biomass <- mask(GEDI_raster, mature_mask)

total_prec <- app(prec, sum)
all_data <- c(mature_biomass, total_prec)

coords <- crds(all_data[[1]], df=FALSE, na.rm=TRUE)
values_stack <- terra::extract(all_data, coords, cells=FALSE, method="simple")
central_df <- values_stack[complete.cases(values_stack), ]

normalized_df <- as.data.frame(lapply(central_df, minMax))


fit <- lm(normalized_df$agbd ~ normalized_df$age)
summary(fit)
plot(normalized_df$agbd ~ normalized_df$age)
abline(fit, col = 'red')

######################################################################
#################        passing into the model     ##################
######################################################################
data <- central_df
# create another list that sets to zero if things are not present

G <- function(pars) {
  # Extract parameters of the model
  E = pars[1]*data$prec + pars[2]*data$temp + pars[3]*data$mat_sum
  LU = pars[8]*data$total_fires * exp(pars[4]*data$ts_fire) + pars[5]*data$pasture + pars[6]*data$other_annual + pars[7]*data$other_perennial
  k = E + LU
  # Prediction of the model
  return ( 1-(exp(-(k))))
}

G <- function(pars) {
  # Extract parameters of the model
  LU = pars[1] * data$total_fires 
  E = pars[2]*data$prec + pars[3]*data$temp + pars[4]*data$mat_sum
  k = E + LU
  # Prediction of the model
  return ( 1-(exp(-(k)))) #normalized_df$Bmax *   ) 
}

G <- function(pars) {
  # Extract parameters of the model
  k = pars[1] * data$age# + pars[2]*data$prec 
  # Prediction of the model
  return ( 1-(exp(-(k)))) #normalized_df$Bmax *   ) 
}

pars = c(0.05,0.05)#, 0.05, 0.05,0.05) #,0.005,0.0005,0.05,0.5)
Gpred <- G(pars)
Gpred

# asymptote is being repeated.
# have so that the different selections have variance.
# normalize the variables so they dont vary wildly. divide the parameter values.

####################
#finding the right parameters

NLL = function(pars) {
  # Values prediced by the model
  if(pars[length(pars)] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  # before returning likelihood value, check if it is NA or inf or something problematic
  Gpred = G(pars)
  #print(Gpred)
  # Negative log-likelihood 
  fun = -sum(dnorm(x = data$agbd - Gpred, mean = 0, sd = pars[length(pars)], log = TRUE), na.rm = TRUE)
  return(fun)
}

o = optim(par = pars, fn = NLL, hessian = FALSE)
print(o)

# meth0 = c("Nelder-Mead", "BFGS", "CG")
# for (i in meth0){
# o = optim(par = pars, fn = NLL, control = list(parscale = abs(pars)), 
#             hessian = FALSE, method = i)
# print(i)
# print(o)
# }

pred = G(o$par[1:length(o$par)])

#jpeg('model2.jpeg')
plot(data$agbd, pred, abline(0,1))
#dev.off()





tst <- pred[pred < 0]
# converting the cells from one to the other is not needed here since I am using two rasters
# add total of b neighboring cells

o$par[1:length(o$par)]
# do general GLM or GAM to see if the environmental variables predict height.
# investigate how the environmental predictors work with mature forests.
# sum expected growth over each year to get final biomass. integration
  # shape of curve will change - rate of growth decreases.
  # think of the asymptote - 

# why are these values negative - check the functional form for what happens after -1.

abline(fit)



# are the likelihoods varying?
# if the surface 
