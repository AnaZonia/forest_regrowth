library(terra)
library(doParallel)

set.seed(0)
mature_biomass <- rast("./data/mature_biomass_1k.tif")
sec_dist <- rast("./data/distance.tif")

sec_dist_coords <- xyFromCell(sec_dist, which(!is.na(values(sec_dist))))
sec_dist_coords <- sec_dist_coords[sample(nrow(sec_dist_coords), 25000, replace = FALSE), ]

distances <- extract(sec_dist, sec_dist_coords)

#  Calculate mean biomass within a buffer
mean_biomass_within_buffer <- function(coord, dist, biomass_raster) {
    x <- coord[1]
    y <- coord[2]
    point <- vect(cbind(x, y), crs = crs(biomass_raster))
    buffer_zone <- buffer(point, width = dist)

    masked_biomass <- mask(biomass_raster, buffer_zone)

    mean_biomass <- global(masked_biomass, fun = mean, na.rm = TRUE)$mean
    return(mean_biomass)
}

# Set up parallel processing
registerDoParallel(cores = 20)

# cl <- makeCluster(num_cores)
# clusterExport(cl, c(
#     "mean_biomass_within_buffer", "sec_dist_coords", "distances",
#     "mature_biomass_matrix", "res_x", "res_y", "extent"
# ))

print(Sys.time())

nearest_biomass <- foreach(i = 1:nrow(sec_dist_coords), .combine = "c") %dopar% {
    print(i)
    mean_biomass_within_buffer(
        sec_dist_coords[i, ], distances[i, ], mature_biomass
    )
}


# Create a new raster to store the results
result <- sec_dist
result[which(!is.na(values(sec_dist)))] <- nearest_biomass
writeRaster(result, "nearest_biomass.tif", overwrite = TRUE)