library(terra)
library(foreach)
library(doParallel)
library(parallel)

set.seed(0)
mature_biomass <- rast("./data/mature_biomass_1k.tif")
sec_dist <- rast("./data/distance.tif")

sec_dist_coords <- xyFromCell(sec_dist, which(!is.na(values(sec_dist))))
sec_dist_coords <- sec_dist_coords[sample(nrow(sec_dist_coords), 30000, replace = FALSE), ]

distances <- extract(sec_dist, sec_dist_coords)

wrapped_biomass <- wrap(mature_biomass, proxy = TRUE)

# Define the function to calculate mean biomass within a buffer
mean_biomass_within_buffer <- function(coord, dist, biomass_raster) {
    x <- coord[1]
    y <- coord[2]
    point <- vect(cbind(x, y), crs = crs(biomass_raster))
    buffer_zone <- buffer(point, width = dist)

    # Mask the biomass raster with the buffer zone
    masked_biomass <- mask(biomass_raster, buffer_zone)

    # Calculate the mean biomass value within the buffer zone
    mean_biomass <- global(masked_biomass, fun = mean, na.rm = TRUE)$mean
    return(mean_biomass)
}

# Function to be used with parLapply
f <- function(i, coord, dist, wrapped_biomass) {
    biomass_raster <- unwrap(wrapped_biomass)
    result <- mean_biomass_within_buffer(coord, dist, biomass_raster)
    return(result)
}

# Set up parallel processing
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c(
    "mean_biomass_within_buffer", "sec_dist_coords",
    "distances", "wrapped_biomass", "f"
))
clusterEvalQ(cl, library(terra))


# Use parLapply for parallel processing
nearest_biomass <- parLapply(cl, 1:nrow(sec_dist_coords), function(i) {
    f(i, sec_dist_coords[i, ], distances[i, ], wrapped_biomass)
})

# Stop the cluster
stopCluster(cl)

# Combine results
nearest_biomass <- unlist(nearest_biomass)

# Create a new raster to store the results
result <- sec_dist
result[which(!is.na(values(sec_dist)))] <- nearest_biomass

# Save the result raster
writeRaster(result, "result.tif", overwrite = TRUE)

# num_cores <- parallel::detectCores() - 1
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)

# # Export necessary objects to the cluster
# clusterExport(cl, c("mean_biomass_within_buffer", "mature_biomass", "sec_dist_coords", "distances"))

# print(Sys.time())

# nearest_biomass <- foreach(i = 1:nrow(distances), .combine = "c", .packages = "terra",
#                             .export = c("r"),
#                             .inorder = TRUE) %dopar% {
#     mean_biomass_within_buffer(
#         sec_dist_coords[i, ],
#         distances[i, ],
#         r
#     )
# }

# print(Sys.time())

# # Shut down parallel processing
# stopCluster(cl)

# Save the results as a raster file
output_raster <- sec_dist
values(output_raster) <- NA
cell_indices <- cellFromXY(sec_dist, sec_dist_coords)
values(output_raster)[cell_indices] <- nearest_biomass
writeRaster(output_raster, "nearest_biomass.tif", overwrite = TRUE)
