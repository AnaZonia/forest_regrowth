

vars <- c("tmax", "tmin")
var <- 'tmax'
for (var in vars){
  list_clim <- mk_list(var)
  # read the files
  raster_clim <- lapply(list_clim, rast)

  yearly <- c()
  for (i in seq(1, length(raster_clim), 12)){ # 12 months in a year, 408 months
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
}

#resolution is about 4.5 km.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  Soil ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  AET - TerraClimate ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
library(terra) # handling spatial data
library(tidyverse)
setwd("/home/aavila/forest_regrowth")
aet_list <- as.list(list.files(path = './water_deficit_terraclimate', pattern = "*.tif", full.names = TRUE))

raster_clim <- lapply(aet_list, rast)

yearly <- c()
for (i in seq(1, length(raster_clim), 12)){ # 12 months in a year, 408 months
  tst <- raster_clim[[i]]
  print(i)
  for (j in (i+1):(i+11)){
    tst <- tst + raster_clim[[j]]
    print(j)
  }
  yearly <- c(yearly, tst)
}

monthly <- c()
for (i in 1:12){ # 12 months in a year, 408 months
  tst <- raster_clim[[i]]
  print(i)
  for (j in seq(i, length(raster_clim), 12)){
    tst <- tst + raster_clim[[j]]
  }
  monthly <- c(monthly, tst)
}


monthly <- rast(monthly)
monthly[monthly == 0] <- NA
monthly <- terra::trim(monthly)
plot(monthly)
yearly <- rast(yearly)
yearly[yearly == 0] <- NA
yearly <- terra::trim(yearly)
plot(yearly)
plot(yearly[[19:35]])

total_years <- sum(yearly)

cropped <- crop(yearly, regrowth)
regrowth <- crop(regrowth, cropped)

cropped_resampled <- resample(total_sum, regrowth, method='near')
raster_clim_masked <- mask(cropped_resampled, regrowth)
writeRaster(raster_clim_masked, filename='cwd_totalsum_santoro.tif')
cropped_resampled <- resample(yearly, regrowth, method='near')
raster_clim_masked <- mask(cropped_resampled, regrowth)
writeRaster(raster_clim_masked, filename='aet_yearly_santoro.tif')

cropped <- crop(monthly, regrowth)
regrowth <- crop(regrowth, cropped)
cropped_resampled <- resample(monthly, regrowth, method='near')
raster_clim_masked <- mask(cropped_resampled, regrowth)
writeRaster(raster_clim_masked, filename='cwd_monthly_santoro.tif')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  CWD - Chave ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cwd <- rast('/home/aavila/forest_regrowth/heinrich_poorter/CWD_poorter.tif')
cwd <- resample(cwd, regrowth)
rasters <- list(regrowth_paper, santoro_raster, cwd)
fire <- rast('./model_ready_rasters/0000000000-0000095232_fire_history_santoro.tif')
last_LU <- rast('last_LU_0000000000-0000095232_santoro.tif')
last_LU <- crop(last_LU, all_rasters)
fire <- crop(fire, all_rasters)

all_rasters <- c(all_rasters, last_LU, fire)

rast_to_df <- function(raster){
  coords <- crds(raster, df=FALSE, na.rm=TRUE)
  values_stack <- terra::extract(raster, coords, cells=FALSE, method="simple")
  central_df <- values_stack[complete.cases(values_stack), ]
  return(central_df)
}

data <- rast_to_df(all_rasters)
colnames(data) <- c('age', 'agbd', 'cwd', 'last_LU', 'fire_total', 'last_fire')

data2 <- data[data$agbd>0,]
saveRDS(data, 'santoro_cwd_fire_LU.rds')

data <- lapply(lu_tile, rast_to_df)
data_raw <- bind_rows(data)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########  Mature forest ##########
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# this is a long term average.

# mean biomass of surrounding mature forests around 3km (~100 pixels of 30m)
biomass_range <- 101 # window size
mature_biomass <- mask(santoro_raster, mature_mask)  # get only biomass of mature forests
# mature_total_biomass <- focal(mature_biomass, biomass_range, sum, na.rm = TRUE)

regrowth_mat <- as.matrix(regrowth)

coordinates <- which(!is.na(regrowth_mat), TRUE)
center_coordinates_reg <- coordinates[,1]

r <- rast(ncols=36, nrows=18)
r[500] <- 1
r <- as.matrix(r, wide = TRUE)

coordinates <- which(!is.na(r), TRUE)

# Define the circular radius (distance) and buffer size
radius <- 10
buffer_size <- 2 * radius + 1  # Buffer size including the central cell

n_rows <- nrow(r)
n_cols <- ncol(r)
row_indices <- rep(1:n_rows, each = n_cols)
col_indices <- rep(1:n_cols, times = n_rows)
cell_coordinates <- cbind(row_indices, col_indices)
distances <- sqrt((cell_coordinates[, 1] - coordinates[1])^2 + (cell_coordinates[, 2] - coordinates[2])^2)
distance_matrix <- matrix(distances, ncol = n_cols, byrow = TRUE)

R <- matrix(1:100, nrow = 20)
row <- 14
col <- 32
radius <- 3

x <- seq(1, nrow(R))
y <- seq(1, ncol(R))
xx <- rep(x, each = length(y))
yy <- rep(y, length(x))
distances <- sqrt((xx - row)^2 + (yy - col)^2)
indices <- which(distances <= radius)
result <- r[within_radius_mask]

within_radius_mask <- which(distance_matrix <= radius)

mask <- matrix(0, nrow = nrow(r), ncol = ncol(r))
mask[within_radius_mask] <- 1
mask <- rast(mask)
plot(mask)

biom <- matrix(1:(nrow(r)*ncol(r)), nrow = nrow(r), ncol = ncol(r))

values <- biom[within_radius_mask]

result <- matrix(NA, nrow = nrow(r), ncol = ncol(r))
result[coordinates[1], coordinates[2]] <- mean(values)

# Load required libraries
library(terra)
library(parallel)

tiles <- makeTiles(trim(mature_biomass), c(200,200), na.rm = TRUE, overwrite=TRUE)  # Adjust nx and ny to desired block size
tiles <- lapply(tiles, rast)

# function intakes buffer_size in meters and mean_biomass_in_x_pixels is self explanatory
calc_buffer_and_sum <- function(tile, buffer_size, mean_biomass_in_x_pixels) {
  mature_total_biomass <- focal(mature_biomass_msk, mean_biomass_in_x_pixels, mean, na.rm = TRUE)
  mature_biomass_msk <- mask(mature_total_biomass, rast(tile))

}

buffer_mask <- as.mask(dist)
mask_raster <- rast(r, mask = dist)

library(parallel)
tiles_processed <- mclapply(tiles, focal, 51, mean, na.rm = TRUE, ncores = 20)
saveRDS(tiles_processed, file = ".rds")

library(terra)

r <- rast(ncols=36, nrows=18)
r[500] <- 1
b <- buffer(r, width=5000000) 
plot(b)

s <- rast(ncols=36, nrows=18)
values(s) <- runif(ncell(s))

tst <- mask(s, b)


plot(tst)
