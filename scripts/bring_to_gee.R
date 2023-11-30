

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



#---------------------
# if age is 1 in 2018, it was not forest in 2017
# which means I want category in 2017
# if age is 33 in 2018, it was not forest in 1985
# which means I want category in 1985

get_last_layer_val <- function(val){
  mask_val <- val_toNA(regrowth_mask, val, TRUE)
  year_val <- mask(tmp_rasters[[val]], mask_val)
  return(year_val)
}

last_LU_mrgd <- get_last_layer_val(1)
for (val in 2:length(tmp_rasters)){
  print(val)
  last_LU_mrgd <- merge(last_LU_mrgd, get_last_layer_val(val))
}

writeRaster(last_LU_mrgd, 'last_LU.tif')

#################################################################################

########## LAND USE ##########
# VARIABLES OBTAINED
# number of years under each land use tnrype
# time since last observation of each land use type

# INDEX ## 3 = forest
# 15 = pasture
# 39 = soy
# 46 = coffee
# 20 = sugar cane
# 41 = other annual crop
# 48 = other perennial crop

# total years under each land use type
calc_total_yrs <- function(masked_brick, val){
  masked_brick[masked_brick != val] <- 0
  total_past_years <- sum(masked_brick)
  return(total_past_years/val)
}

pasture_total <- calc_total_yrs(lulc_brick_masked, 15)
soy_total <- calc_total_yrs(lulc_brick_masked, 39)
coffee_total <- calc_total_yrs(lulc_brick_masked, 46)
sugar_total <- calc_total_yrs(lulc_brick_masked, 20)
other_annual_total <- calc_total_yrs(lulc_brick_masked, 41)
other_perennial_total <- calc_total_yrs(lulc_brick_masked, 48)

# not necessary to create a function for last observed;
# just use regrowth mask and ages, and look for the land use in that position.
# we reach the column by tst[[lyr_num]][age]

calc_time_since_lu <- function(masked_brick, val){
  masked_brick_flipped <- masked_brick [[ c(rev(order(names(masked_brick)))) ]]
  lu_instances <- which.lyr(masked_brick == val) # position of layer with last observation of land use type designated by value val
  return(nlyr(masked_brick_flipped) - lu_instances)
}

ts_pasture <- calc_time_since_lu(lulc_brick_masked, 15)
ts_soy <- calc_time_since_lu(lulc_brick_masked, 39)
ts_coffee <- calc_time_since_lu(lulc_brick_masked, 46)
ts_sugar <- calc_time_since_lu(lulc_brick_masked, 20)
ts_other_annual <- calc_time_since_lu(lulc_brick_masked, 41)
ts_other_perennial <- calc_time_since_lu(lulc_brick_masked, 48)

# note that there are land use types in last_LU that are not the 4 most common (soy, pasture, other perennial, other annual)
# decided to keep them for now but can be removed later with this:
# last_LU[last_LU != 15 & last_LU != 39 & last_LU != 41 & last_LU != 48] <- NA

lulc <- c(last_LU, ts_pasture, ts_soy, ts_other_annual, ts_other_perennial,
pasture_total, soy_total, other_annual_total, other_perennial_total)
names(lulc) <- c('last_LU', 'ts_pasture', 'ts_soy', 'ts_other_annual', 'ts_other_perennial',
'pasture_total', 'soy_total', 'other_annual_total', 'other_perennial_total')

writeRaster(lulc, paste0('./model_ready_rasters/', location, '_lulc_history_santoro.tif'))



# --------
