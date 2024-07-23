library(terra)
library(parallel)

set.seed(0)
mature_biomass <- rast("./data/mature_biomass_1k.tif")
sec_dist <- rast("./data/distance.tif")

# Convert rasters to matrices
mature_biomass_matrix <- as.matrix(mature_biomass)
sec_dist_matrix <- as.matrix(sec_dist)

# Get raster properties
res_x <- xres(mature_biomass)
res_y <- yres(mature_biomass)
extent <- ext(mature_biomass)

sec_dist_coords <- xyFromCell(sec_dist, which(!is.na(values(sec_dist))))
sec_dist_coords <- sec_dist_coords[sample(nrow(sec_dist_coords), 30000, replace = FALSE), ]

distances <- extract(sec_dist, sec_dist_coords)

# Define the function to calculate mean biomass within a buffer
mean_biomass_within_buffer <- function(coord, dist, biomass_matrix, res_x, res_y, extent) {
    x <- coord[1]
    y <- coord[2]

    # Calculate buffer in matrix indices
    buffer_cells <- ceiling(dist / res_x)

    # Calculate matrix indices for the buffer area
    row <- nrow(biomass_matrix) - floor((y - extent[3]) / res_y)
    col <- floor((x - extent[1]) / res_x) + 1

    row_min <- max(1, row - buffer_cells)
    row_max <- min(nrow(biomass_matrix), row + buffer_cells)
    col_min <- max(1, col - buffer_cells)
    col_max <- min(ncol(biomass_matrix), col + buffer_cells)

    # Extract the buffer area
    buffer_area <- biomass_matrix[row_min:row_max, col_min:col_max]

    # Calculate mean biomass
    mean_biomass <- mean(buffer_area, na.rm = TRUE)
    return(mean_biomass)
}

# Set up parallel processing
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c(
    "mean_biomass_within_buffer", "sec_dist_coords", "distances",
    "mature_biomass_matrix", "res_x", "res_y", "extent"
))

# Use parLapply for parallel processing
nearest_biomass <- parLapply(cl, 1:nrow(sec_dist_coords), function(i) {
    mean_biomass_within_buffer(sec_dist_coords[i, ], distances[i, ], mature_biomass_matrix, res_x, res_y, extent)
})

# Stop the cluster
stopCluster(cl)

# Combine results
nearest_biomass <- unlist(nearest_biomass)

# Create a new raster to store the results
result <- sec_dist
result[which(!is.na(values(sec_dist)))] <- nearest_biomass

# Save the result raster
writeRaster(result, "path/to/result.tif", overwrite = TRUE)
