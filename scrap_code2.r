
### mature forest confusions
setwd("/home/aavila/forest_regrowth")

tst <- rast('./heinrich_poorter/CWD.tif/CWD.tif')


  yearly <- c()
  for (i in seq(1, length(raster_clim), 12)){ # 12 months in a year, 408 months
    tst <- raster_clim[[i]]
    print(i)
    for (j in (i+1):(i+11)){
      tst <- tst + raster_clim[[2]]
      print(j)
    }
    yearly <- c(yearly, tst)
  }

# this is a long term average.

# mean biomass of surrounding mature forests around 3km (~100 pixels of 30m)
 biomass_range <- 101 # window size
 mature_biomass <- mask(santoro_raster, mature_mask)  # get only biomass of mature forests
# mature_total_biomass <- focal(mature_biomass, biomass_range, sum, na.rm = TRUE)

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
#> terra 1.7.39
par(mfrow = c(2, 1))

r <- rast(ncols=36, nrows=18)
r[500] <- 1
b <- buffer(r, width=5000000) 
plot(b)

s <- rast(ncols=36, nrows=18)
values(s) <- runif(ncell(s))

# set FALSE values to NA
# b[!b] <- NA
# or change maskvalues arg to match the mask raster:
tst <- mask(s, b, maskvalues = FALSE)
plot(tst)




data <- readRDS('santoro_ESA_alldata.rds')

data <- data[data$last_LU %in% c(15, 41, 48),]
data$last_LU <- factor(data$last_LU)
dummy_LU <- as.data.frame(model.matrix(~ data$last_LU - 1))
names(dummy_LU) <- c('pasture', 'other_annual', 'other_perennial')
data <- data[,-7]

minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
data <- as.data.frame(lapply(data, minMax))

#data <- as.data.frame(scale(data))

data <- cbind(data, dummy_LU)
