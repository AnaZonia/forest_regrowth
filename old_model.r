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

#tmin <- rast('./model_ready_rasters/tmin_0000000000-0000095232.tif')
#tmax <- rast('./model_ready_rasters/tmax_0000000000-0000095232.tif')
#temp <- tmax-tmin
#writeRaster(temp, paste0('./model_ready_rasters/temp_0000000000-0000095232.tif'))

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
