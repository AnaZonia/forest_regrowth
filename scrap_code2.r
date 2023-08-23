
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

